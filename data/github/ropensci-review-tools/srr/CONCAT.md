<!-- badges: start -->

[![R build
status](https://github.com/ropensci-review-tools/srr/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci-review-tools/srr/actions)
[![codecov](https://codecov.io/gh/ropensci-review-tools/srr/branch/main/graph/badge.svg)](https://codecov.io/gh/ropensci-review-tools/srr)
[![Project Status:
Concept](https://www.repostatus.org/badges/latest/concept.svg)](https://www.repostatus.org/#concept)
<!-- badges: end -->

<!-- README.md is generated from README.Rmd. Please edit that file -->

# srr

“srr” stands for **S**oftware **R**eview **R**oclets, and is
[rOpenSci](https://ropensci.org)’s package for extending documentation
to include additional components specific to the software review
process. The package currently facilitates documenting how statistical
software complies with our collections of [Statistical Software
Standards](https://stats-devguide.ropensci.org/standards.html). Before
proceeding, the answer to an important question: **[What is a
“roclet”](https://github.com/r-lib/roxygen2/issues/1086)?**

-   A roclet is an object used by the
    [`roxygen2`](https://roxygen2.r-lib.org) package to convert
    [`roxygen2`](https://roxygen2.r-lib.org)-style documentation lines
    into some desired form of output.

## Why then?

This package currently serves to aid developers and reviewers of
statistical software in aligning their software against our extensive
[lists of
standards](https://stats-devguide.ropensci.org/standards.html). In
acknowledgement of [Colin Gillespie](https://github.com/csgillespie)’s
sentiments expressed in his keynote speech at the [European R Users
Meeting 2020](https://2020.erum.io/program/keynotes-invited-speakers/):

> Standards are good<br> Standards should be strict<br> No-one reads
> standards

the `srr` package aims to integrate the task of aligning software with
standards within the practice of coding itself, and to make standards
adherence as painless as possible.

## How?

The [`roxygen2`](https://roxygen2.r-lib.org) package parses all
documentation lines from all files in the `R/` directory of a package
which begin with `#'`. Special tags beginning with `@`, such as `@param`
or `@export`, may follow these symbols, and roclets define what is done
with different kinds of tags. The
[`roxygen2`](https://roxygen2.r-lib.org) package includes roclets to
process a number of tags; the `srr` package implements custom roclets to
process several additional tags for use with
[rOpenSci](https://ropensci.org)’s software review systems.

At present, the package only contains roclets and associated functions
to help those developing and reviewing packages submitted to rOpenSci’s
system for [Statistical Software
Review](https://stats-devguide.ropensci.org/). The functions are mostly
intended to ease alignment and assessment of software against the
standards detailed in the [main project
book](https://stats-devguide.ropensci.org/standards.html) (from here on
referred to as the “SSR Book”).

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
install.packages("srr")
```

Alternatively, the package can be installed by running one of the
following lines:

``` r
remotes::install_github ("ropensci-review-tools/srr")
pak::pkg_install ("ropensci-review-tools/srr")
```

and loaded for use with,

``` r
library (srr)
```

## Overview

Both this `README`, and the main package vignette, describe the
functionality of the package in the specific context of the statistical
software review project. Both the roclet and all functions intended for
use in this context are prefixed with `srr_stats_`. The remainder of
this document is in two main sections. If you’re developing a statistics
package for submission to our peer review system, keep straight on
reading. If you’ve been invited to review a package, you may skip the
following section and just read the subsequent section. The general
procedures for both developers and reviewers are described at length in
the [SSR book](https://stats-devguide.ropensci.org/standards.html), with
this `README` intended to provide supporting technical details.

Note that the `srr` package can be applied only within the working
directory of a package. There are no `package` or `path` arguments to
allow functions to be applied to packages anywhere other than in the
current working directory.

## For Package Developers

People intending to develop packages for submission to our system for
peer reviewing statistical software will need to follow the following
general steps. Note that, while the `srr` package has a few functions
which developers may call directly to aid their submission process, most
functionality of this package is implemented via custom [`roxygen2`
“roclets”](https://roxygen2.r-lib.org). The third of the following steps
describes how to link your package with `srr` in order to use these
roclets.

1.  Ensure that your package successfully passes all
    [`autotest`](https://github.com/ropensci-review-tools/autotest)
    tests, potentially including setting `test = FALSE` flags to switch
    off any tests you consider not to be applicable to your package. For
    details, please see the [package documentation for
    `autotest`](https://docs.ropensci.org/autotest/).

2.  Decide which of our in-scope categories of statistical software best
    describe your package. The function
    [`srr_stats_categories()`](https://docs.ropensci.org/srr/reference/srr_stats_categories.html)
    provides a list of currently developed categories for which
    standards have been developed, along with links to the online
    standards for each category:

    ``` r
    srr_stats_categories ()$title
    ```

        ## [1] "General"                                                        
        ## [2] "Bayesian"                                                       
        ## [3] "EDA"                                                            
        ## [4] "Machine Learning"                                               
        ## [5] "Regression and Supervised Learning"                             
        ## [6] "Spatial"                                                        
        ## [7] "Time Series"                                                    
        ## [8] "Dimensionality Reduction, Clustering, and Unsupervised Learning"

    That function also returns links to the full descriptions of each
    category in the [main project
    book](https://stats-devguide.ropensci.org/standards.html). Any
    software within one or more of these categories may be considered
    for review.

3.  Enable your package to use the `srr_stats` roclets by modifying the
    package’s `DESCRIPTION` file so that the `Roxygen` line looks like
    this:

    ``` r
    Roxygen: list(markdown = TRUE, roclets = c ("namespace", "rd", "srr::srr_stats_roclet"))
    ```

    That will load the [“roclet”](https://roxygen2.r-lib.org) used by
    this package to process the documentation of statistical standards
    within your actual code. Note that you do not need to add, import,
    or depend upon the `srr` package anywhere else within the
    `DESCRIPTION` file or your actual package.

4.  Load the `srr` package and generate lists of standards within your
    package’s `/R` folder by running,
    [`srr_stats_roxygen(category = c("<my-category-1>", "<my-category-2>"))`](https://docs.ropensci.org/srr/reference/srr_stats_roxygen.html).
    This will by default create a new file called by default
    `R/srr_stats_standards.R`, the first few lines of which will look
    like this:

        ## [1] "#' srr_stats"                                                                 
        ## [2] "#'"                                                                           
        ## [3] "#' All of the following standards initially have `@srrstatsTODO` tags."       
        ## [4] "#' These may be moved at any time to any other locations in your code."       
        ## [5] "#' Once addressed, please modify the tag from `@srrstatsTODO` to `@srrstats`,"
        ## [6] "#' or `@srrstatsNA`, ensuring that references to every one of the following"

    The file will contain a list of all standards from your nominated
    categories. This file may be renamed, and the individual items moved
    to other locations in other files, but all nominated standards
    should remain in [`roxygen2`](https://roxygen2.r-lib.org) blocks
    somewhere in your source code.

    The `@srrstatsVerbose` line defines a variable which may be used to
    suppress output from the `srrstats` roclet when updating package
    documentation (by setting to `FALSE`). After that comes the list of
    standards, each of which is prefixed by a
    [`roxygen2`](https://roxygen2.r-lib.org) tag, `@srrstatsTODO`. A
    package can only be submitted once all of these `TODO` items have
    been addressed via one of the options described in the following two
    items.

5.  A standard may be addressed by moving the item in the
    `srr-stats-standards.R` file (or wherever you’ve chosen to list
    these within your own package) to one or more places in your code
    where these standards have been addressed. In doing so, the
    [`roxygen2`](https://roxygen2.r-lib.org) tag should be changed from
    `@srrstatsTODO` to `@srrstats`, and the text which initially lists
    the actual standard should be changed to provide a brief description
    of how that standard has been met. Tags for one particular standard
    may be repeated in multiple places within your code, and we
    encourage locating multiple `@srrstats` tags which refer to a
    particular standard at all locations which directly address that
    standard.

6.  Alternatively, any standards which you consider not applicable to
    your software may remain listed in the templated section of the main
    `srr-stats-standards.R` document (or any alternative location), with
    their tag changed from `@srrstatsTODO` to `@srrstatsNA`, and the
    description of the standard removed and replaced by an explanation
    of why you consider that standard not to be applicable to your
    software. These `@srrstatsNA` tags should be collected together
    within a single `roxygen2` block with a title of `NA_standards`, as
    provided in the initial template generated by the
    [`srr_stats_roxygen()`](https://docs.ropensci.org/srr/reference/srr_stats_roxygen.html)
    function. Any non-applicable standards can then just be moved into
    this block, with their `@srrstatsTODO` tags changed to `@srrstatsNA`

7.  Each time you run
    [`devtools::document()`](https://devtools.r-lib.org/reference/document.html)
    or the equivalent
    [`roxygen2::roxygenise()`](https://roxygen2.r-lib.org/reference/roxygenize.html),
    the roclet will scan your package’s documentation for the state of
    standards, and will generate a summary of the result on your screen.

To help developers understand how to use these roclets, this package
includes a function,
[`srr_stats_pkg_skeleton()`](https://docs.ropensci.org/srr/reference/srr_stats_pkg_skeleton.html),
which will generate a skeleton of a package with several
[`roxygen2`](https://roxygen2.r-lib.org) tags inserted throughout the
code. This function returns the directory where the skeleton package has
been created, so running the following two lines will illustrate the
roclets in action:

``` r
d <- srr_stats_pkg_skeleton ()
roxygen2::roxygenise (d)
```

Note that the skeleton package also includes C++ code in a `src/`
directory, so will be compiled the first time your run
[`roxygensise()`](https://roxygen2.r-lib.org/reference/roxygenize.html).
Running a second time will generate cleaner output from the `srr_stats`
roclets only. The tags included in
[`roxygen2`](https://roxygen2.r-lib.org/) blocks in this skeleton
package may be modified, moved, copied, and changed in any way you like
to help you understand how the roclets work. Simply play around with the
[`roxygen2`](https://roxygen2.r-lib.org/) lines and run
[`roxygensise()`](https://roxygen2.r-lib.org/reference/roxygenize.html)
each time to see the effect. Individual standards may be moved to, and
addressed in, any location including the directories `R/`, `src/`, or
`tests/`, and well as in `.Rmd` documentation files such as `README.Rmd`
or package vignettes. The `srr_stats` roclet associated with this
package is able to parse the various `@srrstats` tags in all of these
locations.

### Places where standards can NOT be inserted

While the `srr` package enables standards compliance to be documented
through inserting `@srrstats` tags in as many locations as possible, in
order to ensure compliance is documented as close as possible to the
point within the code where each standard is addressed, it is not
possible to insert `roxygen2` tags in every type of file. In general,
standards may be inserted in any `.R` or `.Rmd` file, and most types of
files in `src` or `inst/include` directories, as long as they are used
with a package able to convert documentation to a corresponding R file
(such as [`Rcpp`](http://www.rcpp.org/)’s generation of `RcppExports.R`
files which include the C++ documentation).

Tags may generally not be placed in any other kinds of files, including
`.md` files such as `CONTRIBUTING.md`, or other files without extensions
such as `DESCRIPTION`, `NAMESPACE`, or `NEWS`. Standards which are best
addressed in such files must be placed in some other generic location
(such as `R/srr-standards.R`), with a cross-reference to the file in
which they are actually addressed.

## Code of Conduct

Please note that this package is released with a [Contributor Code of
Conduct](https://ropensci.org/code-of-conduct/). By contributing to this
project, you agree to abide by its terms.
# Contributing to `srr`

<!-- This CONTRIBUTING.md is adapted from https://gist.github.com/peterdesmet/e90a1b0dc17af6c12daf6e8b2f044e7c -->

First of all, thanks for considering contributing to `srr`! 👍 It's
people like you that make it rewarding for us - the project maintainers - to
work on `srr`. 😊

`srr` is an open source project, maintained by people who care.

- repo: https://github.com/ropensci-review-tools/srr
- issues: https://github.com/ropensci-review-tools/srr/issues
- new_issue: https://github.com/ropensci-review-tools/srr/issues/new
- website: https://docs.ropensci.org/srr/
- citation: https://ropensci-review-tools.github.io/srr/authors.html
- email: mailto:mark@ropensci.org

## Code of conduct

Please note that this package is released with a [Contributor Code of
Conduct](https://ropensci.org/code-of-conduct/). By contributing to this
project, you agree to abide by its terms.

## How you can contribute

There are several ways you can contribute to this project. If you want to know
more about why and how to contribute to open source projects like this one, see
this [Open Source Guide](https://opensource.guide/how-to-contribute/).

### Share the love ❤️

Think `srr` is useful? Let others discover it, by telling them in person,
via Twitter or a blog post.

Using `srr` for a paper you are writing? Consider [citing
it](https://github.com/ropensci-review-tools/srr/blob/main/inst/CITATION).

### Ask a question ⁉️

Using `srr` and got stuck? [Browse the
documentation](https://docs.ropensci.org/srr/) to see if you can find a
solution. Still stuck? Post your question as an [issue on
GitHub](https://github.com/ropensci-review-tools/srr/issues). While we
cannot offer user support, we'll try to do our best to address it, as questions
often lead to better documentation or the discovery of bugs.

Want to ask a question in private? Contact the package maintainer by
[email](mailto:mark@ropensci.org).

### Propose an idea 💡

Have an idea for a new `srr` feature? Take a look at [the
documentation](https://docs.ropensci.org/srr/) and [issues
list](https://github.com/ropensci-review-tools/srr/issues) to see if it
isn't included or suggested yet. If not, suggest your idea as an [issue on
GitHub](https://github.com/ropensci-review-tools/srr/issues/new).
While we can't promise to implement your idea, it helps to:

* Explain in detail how it would work.
* Keep the scope as narrow as possible.

See below if you want to contribute code for your idea as well.

### Report a bug 🐛

Using `srr` and discovered a bug? That's annoying! Don't let others have
the same experience and report it as an [issue on
GitHub](https://github.com/ropensci-review-tools/srr/issues/new) so we can
fix it. A good bug report makes it easier for us to do so, so please:

- Use [the `reprex` package](https://reprex.tidyverse.org) to create a
  reproducible example.
- Include the version of `srr` with the following line in your `reprex`
  code:
  ```
  packageVersion("srr")
  ```

### Improve the documentation 📖

Noticed a typo on the website? Think a function could use a better example?
Good documentation makes all the difference, so your help to improve it is very
welcome!

#### The website

[This website](https://docs.ropensci.org/srr/) is generated with
[`pkgdown`](http://pkgdown.r-lib.org/). That means we don't have to write any
html: content is pulled together from documentation in the code, vignettes,
[Markdown](https://guides.github.com/features/mastering-markdown/) files, the
package `DESCRIPTION` and `_pkgdown.yml` settings. If you know your way around
`pkgdown`, you can [propose a file
change](https://help.github.com/articles/editing-files-in-another-user-s-repository/)
to improve documentation. If not, [report an
issue](https://github.com/ropensci-review-tools/srr/issues/new) and we can
point you in the right direction.

#### Function documentation

Functions are described as comments near their code and translated to
documentation using [`roxygen2`](https://klutometis.github.io/roxygen/). If you
want to improve a function description:

1. Go to `R/` directory in the [code
   repository](https://github.com/ropensci-review-tools/srr/tree/main/R).
2. Look for the file with the name of the function.
3. [Propose a file
   change](https://help.github.com/articles/editing-files-in-another-user-s-repository/)
   to update the function documentation in the roxygen comments (starting with
   `#'`).

### Contribute code 📝

Care to fix bugs or implement new functionality for `srr`? Awesome! 👏
Have a look at the [issue
list](https://github.com/ropensci-review-tools/srr/issues) and leave a
comment on the things you want to work on. See also the development guidelines
below.

## Development guidelines

We try to follow the [GitHub
flow](https://guides.github.com/introduction/flow/) for development.

1. Fork [this repo](https://github.com/ropensci-review-tools/srr/) and
   clone it to your computer. To learn more about this process, see [this
   guide](https://guides.github.com/activities/forking/).
2. If you have forked and cloned the project before and it has been a while
   since you worked on it, [pull changes from the original
   repo](https://help.github.com/articles/merging-an-upstream-repository-into-your-fork/)
   to your clone by using `git pull upstream master`.
3. Open the RStudio project file (`.Rproj`).
4. Make your changes:
    * Write your code.
    * Test your code (bonus points for adding unit tests).
    * Document your code (see function documentation above).
    * Check your code with `devtools::check()` and aim for 0 errors and warnings.
5. Commit and push your changes.
6. Submit a [pull
   request](https://guides.github.com/activities/forking/#making-a-pull-request).

## Code style

The `srr` coding style diverges somewhat from [the commonly used tidyverse
style guide](https://style.tidyverse.org/syntax.html#spacing), primarily
through judicious use of
whitespace, which aims to improve code readability. Code references in
`srr` are separated by whitespace, just like words of text. Just like it
is easier to understand "these three words" than "thesethreewords", code is 
formatted like this:

``` r
these <- three (words (x))
```

and not like this:

``` r
these <- three(words(x))
```

The position of brackets is then arbitrary, and we could also write

``` r
these <- three( words (x))
```

`srr` code opts for the former style, with the natural result that one
ends up writing

```r
this <- function ()
```

with a space between `function` and `()`. You can easily (re-)format your code
to accord with this style [by installing the `spaceout`
package](https://github.com/ropensci-review-tools/spaceout) and running:

```r
styler::style_pkg (style = spaceout::spaceout_style)
```
---
title: "srr"
output:
  md_document:
    variant: gfm

  rmarkdown::html_vignette:
    self_contained: no
---

<!-- badges: start -->
[![R build status](https://github.com/ropensci-review-tools/srr/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci-review-tools/srr/actions)
[![codecov](https://codecov.io/gh/ropensci-review-tools/srr/branch/main/graph/badge.svg)](https://codecov.io/gh/ropensci-review-tools/srr)
[![Project Status:
Concept](https://www.repostatus.org/badges/latest/concept.svg)](https://www.repostatus.org/#concept)
<!-- badges: end -->


<!-- README.md is generated from README.Rmd. Please edit that file -->

# srr

"srr" stands for **S**oftware **R**eview **R**oclets, and is
[rOpenSci](https://ropensci.org)'s package for extending documentation to
include additional components specific to the software review process. The
package currently facilitates documenting how statistical software complies
with our collections of [Statistical Software
Standards](https://stats-devguide.ropensci.org/standards.html). Before
proceeding, the answer to an important question: **[What is
a "roclet"](https://github.com/r-lib/roxygen2/issues/1086)?**

- A roclet is an object used by the [`roxygen2`](https://roxygen2.r-lib.org)
  package to convert [`roxygen2`](https://roxygen2.r-lib.org)-style
  documentation lines into some desired form of output.

## Why then?

This package currently serves to aid developers and reviewers of statistical
software in aligning their software against our extensive [lists of
standards](https://stats-devguide.ropensci.org/standards.html).
In acknowledgement of [Colin Gillespie](https://github.com/csgillespie)'s
sentiments expressed in his keynote speech at the [European R Users Meeting
2020](https://2020.erum.io/program/keynotes-invited-speakers/):

> Standards are good<br>
> Standards should be strict<br>
> No-one reads standards

the `srr` package aims to integrate the task of aligning software with
standards within the practice of coding itself, and to make standards adherence
as painless as possible. 

## How?

The [`roxygen2`](https://roxygen2.r-lib.org) package parses all documentation
lines from all files in the `R/` directory of a package which begin with `#'`.
Special tags beginning with `@`, such as `@param` or `@export`, may follow
these symbols, and roclets define what is done with different kinds of tags.
The [`roxygen2`](https://roxygen2.r-lib.org) package includes roclets to
process a number of tags; the `srr` package implements custom roclets to
process several additional tags for use with [rOpenSci](https://ropensci.org)'s
software review systems.

At present, the package only contains roclets and associated functions to help
those developing and reviewing packages submitted to rOpenSci's system for
[Statistical Software
Review](https://stats-devguide.ropensci.org/).
The functions are mostly intended to ease alignment and assessment of software
against the standards detailed in the [main project
book](https://stats-devguide.ropensci.org/standards.html)
(from here on referred to as the "SSR Book").

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
install.packages("srr")
```

Alternatively, the package can be installed by running one of the following
lines:


```{r remotes, eval = FALSE}
remotes::install_github ("ropensci-review-tools/srr")
pak::pkg_install ("ropensci-review-tools/srr")
```

and loaded for use with,
```{r lib-fakey, eval = FALSE}
library (srr)
```
```{r lib, echo = FALSE, message = FALSE}
devtools::load_all (".", export_all = FALSE)
```

## Overview

Both this `README`, and the main package vignette, describe the functionality
of the package in the specific context of the statistical software review
project. Both the roclet and all functions intended for use in this context are
prefixed with `srr_stats_`. The remainder of this document is in two main
sections. If you're developing a statistics package for submission to our peer
review system, keep straight on reading. If you've been invited to review a
package, you may skip the following section and just read the subsequent
section. The general procedures for both developers and reviewers are described
at length in the [SSR
book](https://stats-devguide.ropensci.org/standards.html),
with this `README` intended to provide supporting technical details.

Note that the `srr` package can be applied only within the working directory of
a package. There are no `package` or `path` arguments to allow functions to be
applied to packages anywhere other than in the current working directory.

## For Package Developers

People intending to develop packages for submission to our system for peer
reviewing statistical software will need to follow the following general steps.
Note that, while the `srr` package has a few functions which developers may
call directly to aid their submission process, most functionality of this
package is implemented via custom [`roxygen2`
"roclets"](https://roxygen2.r-lib.org). The third of the following steps
describes how to link your package with `srr` in order to use these roclets.

1. Ensure that your package successfully passes all
   [`autotest`](https://github.com/ropensci-review-tools/autotest) tests, potentially
   including setting `test = FALSE` flags to switch off any tests you consider
   not to be applicable to your package. For details, please see the [package
   documentation for `autotest`](https://docs.ropensci.org/autotest/).
2. Decide which of our in-scope categories of statistical software best
   describe your package. The function
   [`srr_stats_categories()`](https://docs.ropensci.org/srr/reference/srr_stats_categories.html)
   provides a list of currently developed categories for which standards
   have been developed, along with links to the online standards for each
   category:

    ```{r available}
    srr_stats_categories ()$title
    ```
    That function also returns links to the full descriptions of each category
    in the [main project
    book](https://stats-devguide.ropensci.org/standards.html).
    Any software within one or more of these categories may be considered for
    review. 
3. Enable your package to use the `srr_stats` roclets by modifying the
   package's `DESCRIPTION` file so that the `Roxygen` line looks like this:
    ```{r roxygen, eval = FALSE}
    Roxygen: list(markdown = TRUE, roclets = c ("namespace", "rd", "srr::srr_stats_roclet"))
    ```
    That will load the ["roclet"](https://roxygen2.r-lib.org) used by this
    package to process the documentation of statistical standards within your
    actual code. Note that you do not need to add, import, or depend upon the
    `srr` package anywhere else within the `DESCRIPTION` file or your actual
    package.
4. Load the `srr` package and generate lists of standards within your
   package's `/R` folder by running,
   [`srr_stats_roxygen(category = c("<my-category-1>", "<my-category-2>"))`](https://docs.ropensci.org/srr/reference/srr_stats_roxygen.html).
   This will by default create a new file called by default
   `R/srr_stats_standards.R`, the first few lines of which will look like this:
    ```{r srr_stats_roxygen, echo = FALSE, message = FALSE}
    f <- "./R/srr-stats-standards.R"
    if (!file.exists (f))
        srr_stats_roxygen ()
    rso <- readLines (f)
    head (rso)
    ```
    The file will contain a list of all standards from your nominated
    categories. This file may be renamed, and the individual items moved to
    other locations in other files, but all nominated standards should remain
    in [`roxygen2`](https://roxygen2.r-lib.org) blocks somewhere in your source
    code.

    The `@srrstatsVerbose` line defines a variable which may be used to
    suppress output from the `srrstats` roclet when updating package
    documentation (by setting to `FALSE`). After that comes the list of
    standards, each of which is prefixed by
    a [`roxygen2`](https://roxygen2.r-lib.org) tag, `@srrstatsTODO`. A package
    can only be submitted once all of these `TODO` items have been addressed
    via one of the options described in the following two items.
5. A standard may be addressed by moving the item in the
   `srr-stats-standards.R` file (or wherever you've chosen to list these within
   your own package) to one or more places in your code where these standards
   have been addressed. In doing so, the
   [`roxygen2`](https://roxygen2.r-lib.org) tag should be changed from
   `@srrstatsTODO` to `@srrstats`, and the text which initially lists the
   actual standard should be changed to provide a brief description of how that
   standard has been met. Tags for one particular standard may be repeated in
   multiple places within your code, and we encourage locating multiple
   `@srrstats` tags which refer to a particular standard at all locations which
   directly address that standard.
6. Alternatively, any standards which you consider not applicable to your
   software may remain listed in the templated section of the main
   `srr-stats-standards.R` document (or any alternative location), with their
   tag changed from `@srrstatsTODO` to `@srrstatsNA`, and the description of
   the standard removed and replaced by an explanation of why you consider that
   standard not to be applicable to your software. These `@srrstatsNA` tags
   should be collected together within a single `roxygen2` block with a title
   of `NA_standards`, as provided in the initial template generated by the
   [`srr_stats_roxygen()`](https://docs.ropensci.org/srr/reference/srr_stats_roxygen.html)
   function. Any non-applicable standards can then just be moved into this
   block, with their `@srrstatsTODO` tags changed to `@srrstatsNA` 
7. Each time you run 
    [`devtools::document()`](https://devtools.r-lib.org/reference/document.html)
    or the equivalent
    [`roxygen2::roxygenise()`](https://roxygen2.r-lib.org/reference/roxygenize.html),
    the roclet will scan your package's documentation for the state of
    standards, and will generate a summary of the result on your screen.

To help developers understand how to use these roclets, this package includes a
function,
[`srr_stats_pkg_skeleton()`](https://docs.ropensci.org/srr/reference/srr_stats_pkg_skeleton.html),
which will generate a skeleton of a package with several
[`roxygen2`](https://roxygen2.r-lib.org) tags inserted throughout the code.
This function returns the directory where the skeleton package has been
created, so running the following two lines will illustrate the roclets in
action:

```{r demo, eval = FALSE}
d <- srr_stats_pkg_skeleton ()
roxygen2::roxygenise (d)
```

Note that the skeleton package also includes C++ code in a `src/` directory, so
will be compiled the first time your run
[`roxygensise()`](https://roxygen2.r-lib.org/reference/roxygenize.html).
Running a second time will generate cleaner output from the `srr_stats` roclets
only. The tags included in 
[`roxygen2`](https://roxygen2.r-lib.org/) blocks in this skeleton package may
be modified, moved, copied, and changed in any way you like to help you
understand how the roclets work. Simply play around with the 
[`roxygen2`](https://roxygen2.r-lib.org/) lines and run 
[`roxygensise()`](https://roxygen2.r-lib.org/reference/roxygenize.html) each
time to see the effect. Individual standards may be moved to, and addressed in,
any location including the directories `R/`, `src/`, or `tests/`, and well as
in `.Rmd` documentation files such as `README.Rmd` or package vignettes. The
`srr_stats` roclet associated with this package is able to parse the various
`@srrstats` tags in all of these locations.

### Places where standards can NOT be inserted

While the `srr` package enables standards compliance to be documented through
inserting `@srrstats` tags in as many locations as possible, in order to ensure
compliance is documented as close as possible to the point within the code
where each standard is addressed, it is not possible to insert `roxygen2` tags
in every type of file. In general, standards may be inserted in any `.R` or
`.Rmd` file, and most types of files in `src` or `inst/include` directories, as
long as they are used with a package able to convert documentation to
a corresponding R file (such as [`Rcpp`](http://www.rcpp.org/)'s generation of
`RcppExports.R` files which include the C++ documentation).

Tags may generally not be placed in any other kinds of files, including `.md`
files such as `CONTRIBUTING.md`, or other files without extensions such as
`DESCRIPTION`, `NAMESPACE`, or `NEWS`. Standards which are best addressed in
such files must be placed in some other generic location (such as
`R/srr-standards.R`), with a cross-reference to the file in which they are
actually addressed.

## Code of Conduct

Please note that this package is released with a [Contributor Code of
Conduct](https://ropensci.org/code-of-conduct/). By contributing to this
project, you agree to abide by its terms.
---
title: "Preparing Statistical Software with the srr package"
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
  %\VignetteIndexEntry{Preparing Statistical Software with the srr package}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

`srr` stands for **S**oftware **R**eview **R**oclets, and is an R package
containing roclets for general use in helping those developing and reviewing
packages submitted to rOpenSci. At present, the `srr` package only contains
roclets and associated functions to help those developing and reviewing
packages submitted to rOpenSci's system for [Statistical Software
Review](https://stats-devguide.ropensci.org/). This vignette demonstrates how
developers are intended to use this package to document the alignment of their
software with rOpenSci's [standards for statistical
software](https://stats-devguide.ropensci.org/standards.html).

The main functions of the package are constructed as [`roxygen2`
"roclets"](https://roxygen2.r-lib.org), meaning that these functions are called
each time package documentation is updated by calling
[`devtools::document()`](https://devtools.r-lib.org/reference/document.html) or
equivalently
[`roxygen2::roxygenise()`](https://roxygen2.r-lib.org/reference/roxygenize.html).
The former is just a wrapper around the latter, to enable documentation to
updated by a single call from within the [`devtools`
package](https://devtools.r-lib.org). From here on, this stage of updating
documentation and triggering the `srr` roclets will be referred to as "running
[`roxygenise`](https://roxygen2.r-lib.org/reference/roxygenize.html)."



## 1. The package skeleton

The [`srr_stats_pkg_skeleton()`
function](https://docs.ropensci.org/srr/reference/srr_stats_pkg_skeleton.html)
included with this package differs from
[other](https://stat.ethz.ch/R-manual/R-patched/library/utils/html/package.skeleton.html)
[package](https://usethis.r-lib.org/reference/create_package.html)
[skeleton](https://rdrr.io/cran/Rcpp/man/Rcpp.package.skeleton.html) functions.
Rather than providing a skeleton from which you can construct your own package,
the [`srr_stats_pkg_skeleton()`
function](https://docs.ropensci.org/srr/reference/srr_stats_pkg_skeleton.html)
generates a skeleton to help developers understand this package works. This
skeleton of a package is intended to be modified to help you understand which
kinds of modifications are allowed, and which generate errors. The package is
by default called `"demo"`, is constructed in R's
[`tempdir()`](https://stat.ethz.ch/R-manual/R-devel/library/base/html/tempfile.html),
and looks like this:

```{r dir-tree}
library (srr)
d <- srr_stats_pkg_skeleton (pkg_name = "package")
fs::dir_tree (d)
```

The files listed there mostly exist to illustrate how standards can be included
within code. The format of standards can be seen by examining any of those
files. For example, the `test.R` file looks like this:
```{r test-file-dummy, eval = FALSE}
readLines (file.path (d, "R", "test.R"))
```
```{r test-file, echo = FALSE}
print (readLines (file.path (d, "R", "test.R")), width = 20)
```

That file illustrates some of the [`roxygen2`](https://roxygen2.r-lib.org) tags
defined by this package, and described in detail below. These tags are parsed
whenever package documentation is updated with
[`roxygenise()`](https://roxygen2.r-lib.org/reference/roxygenize.html),
which will produce output like the following:

```{r skeleton-output1-compile, eval = TRUE, echo = FALSE, message = FALSE, results = "hide"}
# need to run once to compile /src so compile output does not appear on 2nd time
pkgbuild::compile_dll (d)
```
```{r skeleton-output-dummy, eval = FALSE, echo = TRUE, collapse = TRUE}
roxygen2::roxygenise (d)
```
```{r skeleton-output, eval = TRUE, echo = FALSE, message = TRUE, collapse = TRUE}
roxygen2::roxygenise (d)
```

The "roclets" contained within this package parse any instances of the
package-specified tags described below, and summarise the output by listing all
locations of each kind of tag. Locations are given as file and line numbers,
and where appropriate, the names of associated functions. We recommend that
developers familiarise themselves with the system by simply modifying any
[`roxygen2`](https://roxygen2.r-lib.org) block in any of the files of this
package skeleton, and running
[`roxygenise()`](https://roxygen2.r-lib.org/reference/roxygenize.html)
to see what happens.

## 2. Enabling `srr` roclets for a package

The "roclets" can be enabled for a package by modifying the `DESCRIPTION` file
so that the `Roxygen` line looks like this:

```{r roxygen, eval = FALSE}
Roxygen: list (markdown = TRUE, roclets = c ("namespace", "rd", "srr::srr_stats_roclet"))
```

That will load the ["roclets"](https://roxygen2.r-lib.org) used by this package
to process standards as documented within your actual code. Note that you do
not need to add, import, or depend upon the `srr` package anywhere else within
the `DESCRIPTION` file. See the `DESCRIPTION` file of the package skeleton for
an example.

## 3. `roxygen2` tags

The `srr` packages recognises and uses the following three tags, each of which
may be inserted into [`roxygen2`](https://roxygen2.r-lib.org) documentation
blocks anywhere in your code. The tags are:

1. `@srrstats` tags to indicate standards which have been addressed. These should
   be placed at the locations within your code where specific standards have
   been addressed, and should include a brief textual description of how the
   code at that location addresses the nominated standard.
2. `@srrstatsTODO` tags to indicate standards which have not yet been addressed.
   The [`srr_stats_roxygen()`
   function](https://docs.ropensci.org/srr/reference/srr_stats_roxygen.html)
   described below will place all standards for your nominated categories in an
   initial file in the `/R` directory of your project. Each of these will have
   a tag of `@srrstatsTODO`. As you address each standard, you'll need to move it
   from that initial location to that point in your code where that standard is
   addressed, and change the tag from `@srrstatsTODO` to `@srrstats`. It may help to
   break the initial, commonly very long, single list of standards into smaller
   groups, and to move each of these into the approximate location within your
   code where these are likely to be addressed. For example, all standards
   which pertain to tests can be moved into the `/tests` (or `/tests/testthat`)
   directory.
3. `@srrstatsNA` tags to indicate standards which you deem not applicable to your
   package. These need to be grouped together in a single block with a title of
   `NA_standards`. Such a block is both included in the package skeleton, and
   also in the output of [`srr_stats_roxygen()`
   function](https://docs.ropensci.org/srr/reference/srr_stats_roxygen.html).
   The numbers of any non-applicable standards can then be moved to this block,
   with a note explaining why these standards have been deemed not to apply to
   your package.


The output illustrated above shows how the `srr` package groups these three
tags together, first collecting the output of standards which have been
addressed (via `@srrstats` tags), then showing output of non-applicable standards
(via `@srrstatsNA` tags), and finally standards which are still to be addressed
(via `@srrstatsTODO` tags).

### 3.1 Format of `srr` documentation

The tags shown above from the skeleton package function `test.R` indicates the
expected format of standards within [`roxygen2`](https://roxygen2.r-lib.org)
documentation blocks. The package is designed to error on attempting to run
[`roxygenise()`](https://roxygen2.r-lib.org/reference/roxygenize.html) with any
inappropriately formatted entries, and so should provide informative messages
to help rectify any problems. `srr` documentation must adhere to the following
formatting requirements:

1. The tags must appear immediately after the
   [`roxygen2`](https://roxygen2.r-lib.org)-style comment symbols. This will
   almost always mean `#' @srrstats` (but see details below, under *Locations of
   `srrstats` documentation*).
2. The standards number(s) must be enclosed within curly braces, such as 
   `#' @srrstats {G1.0}`. Multiple standards may be associated with a single tag by
   including them within the same single pair of curly braces, and separating
   each by a comma, such as `#' @srrstats {G1.0, G1.1}`.
3. Explanatory text may be placed anywhere before or after curly braces
   enclosing standards, such that `#' @srrstats some text {G1.0}` is entirely
   equivalent to `#' @srrstats {G1.0} some text`.
4. Only the first pair of curly braces is parsed to determine standards
   numbers; any subsequent curly braces within one expression will be ignored.
   (Where an expression is everything that comes after one
   [`roxygen2`](https://roxygen2.r-lib.org) tag, and extends until the start of
   the next tag.) The following will accordingly only refer to G1.0, and not
   G1.1:
    ```{r}
   #' @srrstats {G1.0} as well as {G1.1}.
    ```
    The appropriate way to refer to multiple tags is to include them in the one set of curly braces, as shown above, or to use separate tags, like this:
    ```{r}
    #' @srrstats {G1.0} as well as
    #' @srrstats {G1.1}
    ```
5. Any standard may occur arbitrarily many times in any file(s) within a
   package. The only requirement is that any one standard should only be
   associated with one single kind of tag; thus should only be `@srrstats` (a
   standard which has been addressed), *or* `@srrstatsNA` (a standard which is not
   applicable), *or* `@srrstatsTODO` (a standard which has not yet been addressed).

In almost all cases, all tags for a package will be generated by the initial
call to
[`srr_stats_roxygen()`](https://docs.ropensci.org/srr/reference/srr_stats_roxygen.html),
and should simply be moved to appropriate locations within a package's code
without modifying the format.

### 3.2 Locations of `srrstats` documentation

`@srrstats` tags and accompanying text can be placed almost anywhere within a
package, especially in any file in the main `/R`, `/tests`, or `src/`
directories. Within the `/R` directory, tags should be placed only in
[`roxygen2`](https://roxygen2.r-lib.org) blocks. These tags, and all associated
text, will be ignored by the roclets used by
[`roxygen2`](https://roxygen2.r-lib.org) to generate package documentation, and
will only appear on screen output like that shown above, and generated when
running [`roxygenise()`](https://roxygen2.r-lib.org/reference/roxygenize.html).
If tags need to be inserted where there is no
[`roxygen2`](https://roxygen2.r-lib.org/reference/roxygenize.html) block, then
a new block will need to be created, the minimal requirement for which is that
it then have an additional `#' @noRd` tag to suppress generation of
documentation (`.Rd`) files. If that block is associated with a function, the
following two lines will suffice:

```{r doc-ex-with-fn}
#' @srrstats G1.0 This standard belongs here
#' @noRd
myfunction <- function (...) {
    # ...
}
```
If that block is not associated with a function, the documentation can be
followed by a `NULL`, like this:
```{r doc-ex-with-null, eval = FALSE}
#' @srrstats G1.0 This standard belongs here
#' @noRd
NULL
```
Unlike [`roxygen2`](https://roxygen2.r-lib.org), which only processes blocks
from within the main `R/` directory, the `srr` package process blocks from
within other directories too. As these blocks will never be passed through to
the main [`roxygenise()`](https://roxygen2.r-lib.org/reference/roxygenize.html)
function, they need neither `#' @noRd` tags, nor `NULL` definitions where none
otherwise exist. An example is in the `tests/` directory of the package
skeleton:
```{r tests-ex-dummy, eval = FALSE}
readLines (file.path (d, "tests", "testthat", "test-a.R"))
```
```{r tests-ex, echo = FALSE}
print (readLines (file.path (d, "tests", "testthat", "test-a.R")),
       width = 20)
```

While `@srrstats` tags can also be placed in the `src/` directory, the package
currently only parses [`doxygen`](https://github.com/doxygen/doxygen)-style
blocks for code written in C or C++. (Note that `srr` is currently further
restricted to C++ code compiled with
[`Rcpp`](https://cran.r-project.org/package=Rcpp), but will soon be adapted to
work with other C++ interfaces such as [`cpp11`](https://cpp11.r-lib.org).)
These blocks are converted by 
[`Rcpp`](https://cran.r-project.org/package=Rcpp) into 
[`roxygen2`](https://roxygen2.r-lib.org) blocks in a file called
`R/RcppExports.R`, and so need to include an additional `@noRd` tag to
(optionally) suppress generation of `.Rd` documentation. The skeleton package
again gives an example:
```{r src-ex-dummy, eval = FALSE}
readLines (file.path (d, "src", "cpptest.cpp"))
```
```{r src-ex, echo = FALSE}
print (readLines (file.path (d, "src", "cpptest.cpp")), width = 20)
```

### 3.3 Documenting standards for documentation

Many standards refer to general package documentation, particularly in a
`README` document, or in package vignettes. All such documents are presumed to
be written in `.Rmd` format, for which `@srrstats` tags must be included within
distinct code chunks. Again, the skeleton package has an example as follows:
```{r readme-ex}
readLines (file.path (d, "README.Rmd"))
```

Those lines illustrate the expected form. `@srrstats` tags should be within
a single block contained within a dedicated code chunk. `@srrstats` chunks within
`.Rmd` files will generally use `echo = FALSE` to prevent them appearing in the
rendered documents. The [`roxygen2`](https://roxygen2.r-lib.org) lines do not
need to be followed by a `NULL` (or any other non-
[`roxygen2`](https://roxygen2.r-lib.org) statement), although if additional
R code is necessary for any reason, you may also need to add `eval = FALSE`.

### 3.4 `@srrstatsNA` tags for non-applicable standards

While `@srrstatsTODO` and `@srrstats` tags may be placed anywhere within a package,
`@srrstatsNA` tags used to denote non-applicable standards must be placed within
a dedicated [`roxygen2`](https://roxygen2.r-lib.org) block with a title of
`NA_standards`. As described above, both the package skeleton and the file
produced by calling
[`srr_stats_roxygen()`](https://docs.ropensci.org/srr/reference/srr_stats_roxygen.html)
include templates for this block. The following illustrates a minimal form:

```{r na-block, eval = FALSE}
#' NA_standards
#'
#' @srrstatsNA {S3.3} is not applicable
#' @noRd
NULL
```

An `NA_standards` block must end with `NULL`, rather than be associated with
a function definition. There can be multiple `NA_standards` blocks in any
location, enabling these standards to be moved to approximate locations where
they might otherwise have been addressed. (For example, non-applicable
standards referring to tests might all be grouped together in a single
`NA_standards` block in the `tests/` directory.)

## 4. The `srr` workflow

The first step for any developer intending to use this package on the way to
a submission to rOpenSci's project for peer-reviewing statistical software is
to generate the package skeleton described above, and to try any and all
conceivable ways to modify locations, formats, and other properties of the
[`roxygen2`](https://roxygen2.r-lib.org) tags defined in this package, in order
to understand how these tags are used to generate the summary results when
running [`roxygenise`](https://roxygen2.r-lib.org/reference/roxygenize.html).
Following that initial familiarisation, a typical workflow will involve the
following general steps:

1. Automatically download and insert lists of general and category-specific
   standards in your package by running
   [`srr_stats_roxygen()`](https://docs.ropensci.org/srr/reference/srr_stats_roxygen.html)
   (in the main package directory). This will by default generate a file in the
   `R/` directory called `srr-stats-standards.R`, although any alternative name
   can also be passed to the function (or the file can be renamed after it has
   initially been created).
2. Change your package's `DESCRIPTION` file to use the `srr` roclets by adding
   `roclets = "srr::srr_stats_roclet"`) to the `Roxygen` line, as demonstrated at
   the outset, as well as in the package skeleton.
3. Run [`roxygenise`](https://roxygen2.r-lib.org/reference/roxygenize.html) to
   confirm that the roclets work on your package. You should initially see only
   a single list of `@srrstatsTODO` standards.
4. We recommend as a first step cutting-and-pasting standards to approximate
   locations within your package's code where you anticipate these standards
   being addressed. Multiple copies of any one standard may exist in multiple
   locations, so you may also repeat standards which you anticipate will be
   addressed in multiple locations. This should reduce a potentially very long
   initial list of standards down to several distinct groups of hopefully more
   manageable size.
5. Begin addressing the standards by:
    - Ensuring your code conforms;
    - Moving each standard to the one or more location(s) where you think your
      code most directly addresses them; 
    - Modifying the `@srrstatsTODO` tag to `@srrstats`
    - Changing the initial text describing the standard itself to a brief
      description of how your code addresses that standard.
6. Standards which you deem not to be applicable to your package should be
   grouped together in a single [`roxygen2`](https://roxygen2.r-lib.org) block
   with the title `NA_standards` (as described above, and as generated by the
   [`srr_stats_roxygen()`](https://docs.ropensci.org/srr/reference/srr_stats_roxygen.html)
   function).
7. Update your documentation as frequently as you like or need, and use the
   output of the roclets to inform and guide the process of converting all
   `@srrstatsTODO` tags into either `@srrstats` or `@srrstatsNA` tags.

Note that we do **not** recommend copying files from the package skeleton into
your own package for you `srr` documentation. The following lines demonstrate
what happens if you insert actual standards into the package skeleton:

```{r dl-stds-fakey, echo = TRUE, eval = FALSE}
srr_stats_roxygen (category = "regression", # for example
                   filename = file.path (d, "R", "srr-stats-standards.R"),
                   overwrite = TRUE)
```
```{r dl-stds, echo = FALSE, eval = TRUE, collapse = TRUE}
srr_stats_roxygen (category = "regression", # for example
                   filename = file.path (d, "R", "srr-stats-standards.R"),
                   overwrite = TRUE)
```
```{r mixed-tags-fakey, echo = TRUE, eval = FALSE}
roxygen2::roxygenise (d)
```
```{r mixed-tags, echo = FALSE, eval = TRUE, error = TRUE, collapse = TRUE}
roxygen2::roxygenise (d)
```
To ensure all standards are first inserted with `@srrstatsTODO` tags, and that
there are no duplicates with other tags, please use only the
[`srr_stats_roxygen()`](https://docs.ropensci.org/srr/reference/srr_stats_roxygen.html)
function.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/roclet.R
\name{srr_stats_roclet}
\alias{srr_stats_roclet}
\title{srr_stats_roclet}
\usage{
srr_stats_roclet()
}
\value{
A \pkg{roxygen2} roclet
}
\description{
Get values of all \code{srrstats} tags in function documentation
}
\details{
Note that this function should never need to be called directly. It only
exists to enable "@srrstats" tags to be parsed from \pkg{roxygen2}
documentation.
}
\examples{
srr_stats_roclet ()
}
\seealso{
Other roxygen: 
\code{\link{srr_stats_roxygen}()}
}
\concept{roxygen}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/report.R
\name{srr_report}
\alias{srr_report}
\title{Generate report from \code{ssr} tags.}
\usage{
srr_report(path = ".", branch = "", view = TRUE)
}
\arguments{
\item{path}{Path to package for which report is to be generated}

\item{branch}{By default a report will be generated from the current branch
as set on the local git repository; this parameter can be used to specify any
alternative branch.}

\item{view}{If \code{TRUE} (default), a html-formatted version of the report is
opened in default system browser. If \code{FALSE}, the return object includes the
name of a \code{html}-rendered version of the report in an attribute named 'file'.}
}
\value{
(invisibly) Markdown-formatted lines used to generate the final html
document.
}
\description{
Generate report from \code{ssr} tags.
}
\examples{
\dontrun{
path <- srr_stats_pkg_skeleton ()
srr_report (path)
}
}
\concept{report}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dl-standards.R
\name{srr_stats_roxygen}
\alias{srr_stats_roxygen}
\title{Insert standards into code in \pkg{roxygen2} format}
\usage{
srr_stats_roxygen(
  category = NULL,
  filename = "srr-stats-standards.R",
  overwrite = FALSE
)
}
\arguments{
\item{category}{One of the names of files given in the directory contents of
\url{https://github.com/ropensci/statistical-software-review-book/tree/main/standards},
each of which is ultimately formatted into a sub-section of the standards.}

\item{filename}{Name of 'R' source file in which to write
\pkg{roxygen2}-formatted lists of standards.}

\item{overwrite}{If \code{FALSE} (default) and \code{filename} already exists, a dialog
will ask whether file should be overwritten.}
}
\value{
Nothing
}
\description{
Obtain rOpenSci standards for statistical software, along with one or more
category-specific standards, as a checklist, convert to project-specific
\pkg{roxygen2} format, and save in nominated file.
}
\examples{
\dontrun{
path <- srr_stats_pkg_skeleton ()
# contains a few standards; insert all with:
f <- file.path (path, "R", "srr-stats-standards.R")
file.exists (f)
length (readLines (f)) # only 14 lines
srr_stats_roxygen (category = "regression",
                   file = f,
                   overwrite = TRUE)
length (readLines (f)) # now much longer
}
}
\seealso{
Other roxygen: 
\code{\link{srr_stats_roclet}()}
}
\concept{roxygen}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dl-standards.R
\name{srr_stats_categories}
\alias{srr_stats_categories}
\title{Get details of current statistical software categories}
\usage{
srr_stats_categories()
}
\value{
A \code{data.frame} with 3 columns of "category" (the categories to be
submitted to \link{srr_stats_checklist}), "title" (the full title), and
"url".
}
\description{
List all currently available categories and associated URLs to full category
descriptions.
}
\examples{
srr_stats_categories ()
}
\seealso{
Other helper: 
\code{\link{srr_stats_checklist_check}()},
\code{\link{srr_stats_checklist}()},
\code{\link{srr_stats_pkg_skeleton}()},
\code{\link{srr_stats_pre_submit}()}
}
\concept{helper}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/srr-package.R
\docType{package}
\name{srr-package}
\alias{srr}
\alias{srr-package}
\title{srr: 'rOpenSci' Review Roclets}
\description{
Companion package to 'rOpenSci' statistical software review project.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://docs.ropensci.org/srr/}
  \item \url{https://github.com/ropensci-review-tools/srr}
  \item Report bugs at \url{https://github.com/ropensci-review-tools/srr/issues}
}

}
\author{
\strong{Maintainer}: Mark Padgham \email{mark@ropensci.org} (\href{https://orcid.org/0000-0003-2172-5265}{ORCID})

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/checklist.R
\name{srr_stats_checklist_check}
\alias{srr_stats_checklist_check}
\title{Check a completed standards checklist}
\usage{
srr_stats_checklist_check(file)
}
\arguments{
\item{file}{Name of local file containing a completed checklist. Must be a
markdown document in \code{.md} format, not \code{.Rmd} or anything else.}
}
\description{
Correct any potential formatting issues in a completed standards checklist
}
\examples{
\dontrun{
f <- tempfile (fileext = ".md")
srr_stats_checklist (category = "regression", filename = f)
chk <- srr_stats_checklist_check (f)
}
}
\seealso{
Other helper: 
\code{\link{srr_stats_categories}()},
\code{\link{srr_stats_checklist}()},
\code{\link{srr_stats_pkg_skeleton}()},
\code{\link{srr_stats_pre_submit}()}
}
\concept{helper}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pkg-skeleton.R
\name{srr_stats_pkg_skeleton}
\alias{srr_stats_pkg_skeleton}
\title{Make skeleton package to test roclet system}
\usage{
srr_stats_pkg_skeleton(base_dir = tempdir(), pkg_name = "demo")
}
\arguments{
\item{base_dir}{The base directory where the package should be constructed.}

\item{pkg_name}{The name of the package. The final location of this package
will be in \code{file.path(base_dir, pkg_name)}.}
}
\value{
The path to the directory holding the newly created package
}
\description{
Make a dummy package skeleton including 'srr' \pkg{roxygen2} tags which can
be used to try out the functionality of this package. Running the example
lines below which activate the 'srr' roclets, and show you what the output
of those roclets looks like. Feel free to examine the effect of modifying any
of the \verb{@srrstats} tags within the code as identified by running those lines.
}
\examples{
d <- srr_stats_pkg_skeleton (pkg_name = "mystatspkg")
# (capture.output of initial compliation messages)
x <- capture.output (roxygen2::roxygenise (d), type = "output")
}
\seealso{
Other helper: 
\code{\link{srr_stats_categories}()},
\code{\link{srr_stats_checklist_check}()},
\code{\link{srr_stats_checklist}()},
\code{\link{srr_stats_pre_submit}()}
}
\concept{helper}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pre-submit.R
\name{srr_stats_pre_submit}
\alias{srr_stats_pre_submit}
\title{Perform pre-submission checks}
\usage{
srr_stats_pre_submit(path = ".", quiet = FALSE)
}
\arguments{
\item{path}{Path to local repository to check}

\item{quiet}{If 'FALSE', display information on status of package on screen.}
}
\value{
(Invisibly) List of any standards missing from code
}
\description{
Check that all standards are present in code, and listed either as
'@srrstats' or '@srrstatsNA'
}
\examples{
d <- srr_stats_pkg_skeleton ()
# The skeleton has 'TODO' standards, and also has only a few from the full
# list expected for the categories specified there.
srr_stats_pre_submit (d)
}
\seealso{
Other helper: 
\code{\link{srr_stats_categories}()},
\code{\link{srr_stats_checklist_check}()},
\code{\link{srr_stats_checklist}()},
\code{\link{srr_stats_pkg_skeleton}()}
}
\concept{helper}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dl-standards.R
\name{srr_stats_checklist}
\alias{srr_stats_checklist}
\title{Download checklists of statistical software standards}
\usage{
srr_stats_checklist(category = NULL, filename = NULL)
}
\arguments{
\item{category}{One of the names of files given in the directory contents of
\url{https://github.com/ropensci/statistical-software-review-book/tree/main/standards},
each of which is ultimately formatted into a sub-section of the standards.}

\item{filename}{Optional name of local file to save markdown-formatted
checklist. A suffix of \code{.md} will be automatically appended.}
}
\value{
A character vector containing a markdown-style checklist of general
standards along with standards for any additional categories.
}
\description{
Obtain rOpenSci standards for statistical software, along with one or more
category-specific standards, as a checklist, and store the result in the
local clipboard ready to paste.
}
\examples{
\dontrun{
x <- srr_stats_checklist (category = "regression")
# or write to specified file:
f <- tempfile (fileext = ".md")
x <- srr_stats_checklist (category = "regression", filename = f)
}
}
\seealso{
Other helper: 
\code{\link{srr_stats_categories}()},
\code{\link{srr_stats_checklist_check}()},
\code{\link{srr_stats_pkg_skeleton}()},
\code{\link{srr_stats_pre_submit}()}
}
\concept{helper}
