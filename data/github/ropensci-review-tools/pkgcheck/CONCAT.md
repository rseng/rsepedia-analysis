# pkgcheck

<!-- badges: start -->

[![R build
status](https://github.com/ropensci-review-tools/pkgcheck/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci-review-tools/pkgcheck/actions?query=workflow%3AR-CMD-check)
[![gitlab
push](https://github.com/ropensci-review-tools/pkgcheck/workflows/push-to-gitlab/badge.svg)](https://github.com/ropensci-review-tools/pkgcheck/actions?query=workflow%3Apush-to-gitlab)
[![codecov](https://codecov.io/gh/ropensci-review-tools/pkgcheck/branch/main/graph/badge.svg)](https://codecov.io/gh/ropensci-review-tools/pkgcheck)
[![Project Status:
Concept](https://www.repostatus.org/badges/latest/concept.svg)](https://www.repostatus.org/#concept)
<!-- badges: end -->

Check whether a package is ready for submission to
[rOpenSci](https://ropensci.org)’s peer review system. The primary
function collates the output of
[`goodpractice`](https://github.com/mangothecat/goodpractice), including
`R CMD check` results, a number of statistics via the [`pkgstats`
package](https://github.com/ropensci-review-tools/pkgstats), and checks
for package structure expected for rOpenSci submissions. The output of
this function immediately indicates whether or not a package is “Ready
to Submit”.

## Installation

The easiest way to install this package is via the [associated
`r-universe`](https://ropensci-review-tools.r-universe.dev/ui#builds).
As shown there, simply enable the universe with

``` r
options (repos = c (
    ropenscireviewtools = "https://ropensci-review-tools.r-universe.dev",
    CRAN = "https://cloud.r-project.org"
))
```

And then install the usual way with,

``` r
install.packages ("pkgcheck")
```

Alternatively, the package can be installed by running one of the
following lines:

``` r
remotes::install_github ("ropensci-review-tools/pkgcheck")
pak::pkg_install ("ropensci-review-tools/pkgcheck")
```

The package can then loaded for use with

``` r
library (pkgcheck)
```

## Setup

The [`pkgstats`
package](https://github.com/ropensci-review-tools/pkgstats) also
requires the system libraries [`ctags`](https://ctags.io) and [GNU
`global`](https://www.gnu.org/software/global/) to be installed.
Procedures to install these libraries on various operating systems are
described in the [`pkgstats`
README](https://docs.ropensci.org/pkgstats). This package also uses the
[GitHub GraphQL API](https://developer.github.com/v4) which requires a
local GitHub token to be stored with an unambiguous name including
`GITHUB`, such as `GITHUB_TOKEN` (recommended), or `GITHUB_PAT` (for
Personal Authorization Token). This can be obtained from GitHub (via
your user settings), and stored using

``` r
Sys.setenv ("GITHUB_TOKEN" = "<my_token>")
```

This can also be set permanently by putting this line in your
`~/.Renviron` file (or creating this if it does not yet exist). Once
`pkgstats` has been successfully installed, the `pkgcheck` package can
then be loaded via a `library` call:

``` r
library (pkgcheck)
```

## Usage

The package primarily has one function, `pkgcheck`, which accepts the
single argument, `path`, specifying the local location of a git
repository to be analysed. The following code generates a reproducible
report by first downloading a local clone of a repository called
[`srr-demo`](https://github.com/mpadge/srr-demo), which contains the
skeleton of an [`srr` (Software Review Roclets)
package](https://github.com/ropensci-review-tools/srr), generated with
the [`srr_stats_pkg_skeleton()`
function](https://docs.ropensci.org/srr/reference/srr_stats_pkg_skeleton.html):

``` r
mydir <- file.path (tempdir (), "srr-demo")
gert::git_clone ("https://github.com/mpadge/srr-demo", path = mydir)
x <- pkgcheck (mydir)
```

That object has default `print` and `summary` methods. The latter can be
used to simply check whether a package is ready for submission:

``` r
summary (x)
## 
## ── demo 0.0.0.9000 ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
## 
## ✔ Package name is available
## ✖ does not have a 'CITATION' file.
## ✖ does not have a 'codemeta.json' file.
## ✖ does not have a 'contributing' file.
## ✔ uses 'roxygen2'.
## ✔ 'DESCRIPTION' has a URL field.
## ✖ 'DESCRIPTION' does not have a BugReports field.
## ✖ Package has at no HTML vignettes
## ✔ All functions have examples.
## ✖ Package has no continuous integration checks.
## ✖ Package coverage is 0% (should be at least 75%).
## ✔ R CMD check found no errors.
## ✔ R CMD check found no warnings.
## 
## ℹ Current status:
## ✖ This package is not ready to be submitted.
## 
```

A package may only be submitted when the summary contains all ticks and
no cross symbols. (These symbols are colour-coded with green ticks and
red crosses when generated in a terminal; GitHub markdown only renders
them in black-and-white.) The object returned from the `pkgcheck`
function is a complex nested list with around a dozen primary
components. Full information can be obtained by simply calling the
default `print` method by typing the object name (`x`).

## What is checked?

### Summary of Check Results

Calling `summary()` on the object returned by the [`pkgcheck()`
function](https://docs.ropensci.org/pkgcheck/reference/pkgcheck.html)
will generate a checklist like that shown above. This checklist will
also be automatically generated when a package is first submitted to
rOpenSci, and is used by the editors to assess whether to process a
submission. Authors must ensure prior to submission that there are no
red crosses in the resultant list. (In the unlikely circumstances that a
package is unable to pass particular checks, explanations should be
given upon submission about why those checks fail, and why review may
proceed in spite of such failures.)

The full list of checks which packages are expected to pass currently
includes:

1.  Package must use
    [`roxygen2`](https://devguide.ropensci.org/building.html#roxygen2-use)
    for documentation.
2.  Package must have a [`contributing.md`
    file](https://devguide.ropensci.org/collaboration.html#contributing-guide).
3.  Package must have a [`CITATION` file in the `inst`
    directory](https://cran.r-project.org/doc/manuals/r-release/R-exts.html#CITATION-files).
4.  Package must have a [`codemeta.json`
    file](https://devguide.ropensci.org/building.html#creating-metadata-for-your-package).
5.  All exported functions must include examples in their documentation.
6.  Left-assign operators must be used consistently throughout all code
    (so either all `=` or all `<-`, but not a mixture of both).
7.  Package `DESCRIPTION` file must have a “URL” field.
8.  Package `DESCRIPTION` file must have a “BugReports” field.
9.  Package name must be available (or package must already be) on CRAN.
10. Package must have continuous integration tests.
11. Package must have test coverage of at least 75%.
12. No function should have a [cyclomatic
    complexity](https://github.com/MangoTheCat/cyclocomp) of 15 or
    greater.
13. `R CMD check` (implemented via the [`rcmdcheck`
    package](https://r-lib.github.io/rcmdcheck/)) must generate no
    warnings or errors.
14. All statistical standards must be documented, as confirmed by the
    [`srr::srr_pre_submit()`
    function](https://docs.ropensci.org/srr/reference/srr_stats_pre_submit.html).

Several of these are by default only shown when they fail; absence from
a resultant checklist may be interpreted to indicate successful checks.

### Detailed Check Results

Full details of check results can be seen by `print`-ing the object
returned by the [`pkgcheck()`
function](https://docs.ropensci.org/pkgcheck/reference/pkgcheck.html)
(or just by typing the name of this object in the console.) This object
is itself a list including the following items:

``` r
names (x)
```

    ## [1] "pkg"          "info"         "checks"       "meta"         "goodpractice"

The first four of these contain information on the package. The
remainder include:

-   `summary` containing summary statistics on package structure (such
    as lines of code, and numbers of files and functions).
-   `git` containing information and statistics on the `git` repository
    associated with a package (if any), including date of first commit,
    total number of commits and of committers.
-   `srr` summarising information on statistical packages generated by
    the [`srr` package](https://github.com/ropensci-review-tools/srr),
    including a link to a locally-generated HTML report on standards
    compliance.
-   `file_list` as a list of binary flags denoting presence of absence
    of all files required for an R package to be submitted to rOpenSci
-   `fns_have_exs` with a list of any functions which do not have
    documented examples
-   `left_assigns` with a logical value indicating whether or not global
    assignment operators (`<<-`) are used; and a tally of the two types
    of left-assignment operators (`<-` and `=`). One one of these should
    be non-zero, reflecting consistent usage of the same type.
-   `pkgstats` containing several statistics of package structure, along
    with associated percentiles assessed against the entire distribution
    of all current CRAN packages.
-   `network_file` with the path to a local HTML file containing a
    [.visjs](https://visjs.org/) representation of the network of
    relationships between objects (such as functions) within a package,
    within and between each computer language used in the package.
-   `gp` containing the output of the [`goodpractice`
    pacakge](http://mangothecat.github.io/goodpractice/), itself
    including results from:
    -   the [`rcmdcheck` pacakge](https://r-lib.github.io/rcmdcheck) for
        running `R CMD check`
    -   the [`covr` package](https://covr.r-lib.org) for assessing code
        coverage
    -   the [`cyclocomp`
        package](https://github.com/MangoTheCat/cyclocomp) for
        quantifying the cyclomatic complexity of package functions
    -   the [`desc` package](https://github.com/r-lib/desc#readme) for
        analysing the structure of `DESCRIPTION` files
    -   the [`lintr` package](https://github.com/jimhester/lintr)
    -   additional checks applied by the [`goodpractice`
        pacakge](http://mangothecat.github.io/goodpractice/)

Note that results from [`lintr`
package](https://github.com/jimhester/lintr) are **not** reported in the
check summary, and that [`lintr`
results](https://github.com/jimhester/lintr) are reported only in the
detailed results, and have no influence on whether a package passes the
summary checks.

### HTML-formatted check results

Printing `pkgcheck` results to screen is nice, but sometimes it’s even
better to have a nicely formatted full-screen representation of check
results. The package includes an additional function,
[`checks_to_markdown()`](https://docs.ropensci.org/pkgcheck/reference/checks_to_markdown.html),
with a parameter, `render`, which can be set to `TRUE` to automatically
render a HTML-formatted representation of the check results, and open it
in a browser. The formatting differs only slightly from the terminal
output, mostly through the components of
[`goodpractice`](http://mangothecat.github.io/goodpractice/) being
divided into distinct headings, with explicit statements in cases where
components pass all checks (the default screen output inherits directly
from that package, and only reports components which *do not* pass all
checks).

This
[`checks_to_markdown()`](https://docs.ropensci.org/pkgcheck/reference/checks_to_markdown.html)
function returns the report in markdown format, suitable for pasting
directly into a GitHub issue, or other markdown-compatible place. (The
[`clipr` package](https://github.com/mdlincoln/clipr) can be used to
copy this directly to your local clipboard with `write_clip(md)`, where
`md` is the output of `checks_to_markdown()`.)

## Caching and running `pkgcheck` in the background

Running the [`pgkcheck`
function](https://docs.ropensci.org/pkgcheck/reference/pkgcheck.html)
can be time-consuming, primarily because the
[`goodpractice`](https://github.com/mangothecat/goodpractice) component
runs both a full `R CMD check`, and calculates code coverage of all
tests. To avoid re-generating these results each time, the package saves
previous reports to a local cache, in a `pkgcheck` subdirectory of the
location determined by

``` r
rappdirs::user_cache_dir ()
```

As explained in the help file for that function, these locations are:

| System   | location                                                                                                 |
|----------|----------------------------------------------------------------------------------------------------------|
| Mac OS X | `~/Library/Caches/pkgcheck`                                                                              |
| Linux    | `~/.cache/pkgcheck`                                                                                      |
| Win XP   | `C:\\Documents and Settings\\<username>\\Local Settings\\Application Data\\<AppAuthor>\\pkgcheck\\Cache` |
| Vista    | `C:\\Users\\<username>\\AppData\\Local\\<AppAuthor>\\pkgcheck\\Cache`                                    |

You may manually erase the contents of this `pkgcheck` subdirectory at
any time at no risk beyond additional time required to re-generate
contents. This default location may also be over-ridden by setting an
environmental variable named `PKGCHECK_CACHE_DIR`. By default checks
presume packages use `git` for version control, with checks updated only
when code is updated via `git commit`. Checks for packages that do not
use `git` are updated when any files are modified.

The first time
[`pkgcheck()`](https://docs.ropensci.org/pkgcheck/reference/pkgcheck.html)
is applied to a package, the checks will be stored in the cache
directory. Calling that function a second time will then load the cached
results, and so enable checks to be returned much faster. For code which
is frequently updated, such as for packages working on the final stages
prior to submission, it may still be necessary to repeatedly call
[`pkgcheck()`](https://docs.ropensci.org/pkgcheck/reference/pkgcheck.html)
after each modification, a step which may still be inconveniently
time-consuming. To facilitate frequent re-checking, the package also has
a [`pkgcheck_bg()`
function](https://docs.ropensci.org/pkgcheck/reference/pkgcheck_bg.html)
which is effectively identical to the main [`pkgcheck()`
function](https://docs.ropensci.org/pkgcheck/reference/pkgcheck.html),
except it runs in the background, enabling you to continue coding while
checks are running.

The [`pkgcheck_bg()`
function](https://docs.ropensci.org/pkgcheck/reference/pkgcheck_bg.html)
returns a handle to the [`callr::r_bg()`
process](https://callr.r-lib.org/reference/r_bg.html) in which the
checks are running. Typing the name of the returned object will
immediately indicate whether the checks are still running, or whether
they have finished. That handle is itself an [`R6`
object](http://r6.r-lib.org/) with a number of methods, notably
including
[`get_result()`](https://callr.r-lib.org/reference/get_result.html)
which can be used to access the checks once the process has finished.
Alternatively, as soon as the background process, the normal
(foreground) [`pkgcheck()`
function](https://docs.ropensci.org/pkgcheck/reference/pkgcheck.html)
may be called to quickly re-load the cached results.

## Prior Work

[The `checklist` package](https://github.com/inbo/checklist) for
“checking packages and R code”.

## Code of Conduct

Please note that this package is released with a [Contributor Code of
Conduct](https://ropensci.org/code-of-conduct/). By contributing to this
project, you agree to abide by its terms.

## Contributors


<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->

All contributions to this project are gratefully acknowledged using the [`allcontributors` package](https://github.com/ropenscilabs/allcontributors) following the [all-contributors](https://allcontributors.org) specification. Contributions of any kind are welcome!

### Code

<table>

<tr>
<td align="center">
<a href="https://github.com/mpadge">
<img src="https://avatars.githubusercontent.com/u/6697851?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci-review-tools/pkgcheck/commits?author=mpadge">mpadge</a>
</td>
<td align="center">
<a href="https://github.com/maelle">
<img src="https://avatars.githubusercontent.com/u/8360597?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci-review-tools/pkgcheck/commits?author=maelle">maelle</a>
</td>
<td align="center">
<a href="https://github.com/assignUser">
<img src="https://avatars.githubusercontent.com/u/16141871?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci-review-tools/pkgcheck/commits?author=assignUser">assignUser</a>
</td>
<td align="center">
<a href="https://github.com/annakrystalli">
<img src="https://avatars.githubusercontent.com/u/5583057?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci-review-tools/pkgcheck/commits?author=annakrystalli">annakrystalli</a>
</td>
<td align="center">
<a href="https://github.com/noamross">
<img src="https://avatars.githubusercontent.com/u/571752?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci-review-tools/pkgcheck/commits?author=noamross">noamross</a>
</td>
</tr>

</table>


### Issue Authors

<table>

<tr>
<td align="center">
<a href="https://github.com/piyalkarum">
<img src="https://avatars.githubusercontent.com/u/48254643?u=370433a2ace6a030f2551575bc08fa53664fbd8f&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci-review-tools/pkgcheck/issues?q=is%3Aissue+author%3Apiyalkarum">piyalkarum</a>
</td>
<td align="center">
<a href="https://github.com/christophsax">
<img src="https://avatars.githubusercontent.com/u/1390827?u=ce6363f6da758d1bb85987d021cacc34a81c8837&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci-review-tools/pkgcheck/issues?q=is%3Aissue+author%3Achristophsax">christophsax</a>
</td>
<td align="center">
<a href="https://github.com/steffilazerte">
<img src="https://avatars.githubusercontent.com/u/14676081?u=0c1a4469750eba5a6f398c35ba34a2dd600a873e&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci-review-tools/pkgcheck/issues?q=is%3Aissue+author%3Asteffilazerte">steffilazerte</a>
</td>
</tr>

</table>


### Issue Contributors

<table>

<tr>
<td align="center">
<a href="https://github.com/dgkf">
<img src="https://avatars.githubusercontent.com/u/18220321?u=bef717254e5b877159fa712e2b8ad6952c816064&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci-review-tools/pkgcheck/issues?q=is%3Aissue+commenter%3Adgkf">dgkf</a>
</td>
<td align="center">
<a href="https://github.com/cboettig">
<img src="https://avatars.githubusercontent.com/u/222586?u=dfbe54d3b4d538dc2a8c276bb5545fdf4684752f&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci-review-tools/pkgcheck/issues?q=is%3Aissue+commenter%3Acboettig">cboettig</a>
</td>
<td align="center">
<a href="https://github.com/jhollist">
<img src="https://avatars.githubusercontent.com/u/5438539?u=815aa29d708acd16180967c0ffaf81fc64a08bf4&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci-review-tools/pkgcheck/issues?q=is%3Aissue+commenter%3Ajhollist">jhollist</a>
</td>
</tr>

</table>

<!-- markdownlint-enable -->
<!-- prettier-ignore-end -->
<!-- ALL-CONTRIBUTORS-LIST:END -->
# Contributing to pkgcheck

## Opening issues

The easiest way to note any behavioural curiosities or to request any new
features is by opening a [github
issue](https://github.com/ropensci-review-tools/pkgcheck/issues).


## Development guidelines

If you'd like to contribute changes to `pkgcheck`, we use [the GitHub
flow](https://guides.github.com/introduction/flow/index.html) for proposing,
submitting, reviewing, and accepting changes. If you haven't done this before,
there's a nice [overview of git](https://r-pkgs.org/git.html), as well
as [best practices for submitting pull requests](http://r-pkgs.org/git.html#pr-make)
in the R packages book by Hadley Wickham and Jenny Bryan.

The `pkgcheck` coding style diverges somewhat from [the commonly used tidyverse style
guide](https://style.tidyverse.org/syntax.html#spacing), primarily through judicious use of
whitespace, which aims to improve code readability. Code references in
`pkgcheck` are separated by whitespace, just like words of text. Just like it
is easier to understand "these three words" than "thesethreewords", code is not
formatted like this:

``` r
these <- three(words(x))
```
rather like this:

``` r
these <- three (words (x))
```

The position of brackets is then arbitrary, and we could also write

``` r
these <- three( words (x))
```

`pkgcheck` code opts for the former style, with the natural result that one
ends up writing

```r
this <- function ()
```

with a space between `function` and `()`. That's it.

## Adding new checks

New checks are a welcome contribution to `pkgcheck`, for which there is a
[dedicated
vignette](https://docs.ropensci.org/pkgcheck/articles/extending-checks.html).
Please discuss any proposed new checks by opening an issue on the GitHub
repository.


## Code of Conduct

We want to encourage a warm, welcoming, and safe environment for contributing to
this project. See the [code of
conduct](https://ropensci.org/code-of-conduct/) for
more information.
# use_github_action_pkgcheck [plain]

    `file_name` must be a character argument

---

    `file_name` must be a single value

---

    `inputs` must be a named list!

# use_github_action_pkgcheck [ansi]

    [1m[22m[30m[47m`file_name`[49m[39m must be a character argument

---

    [1m[22m[30m[47m`file_name`[49m[39m must be a single value

---

    [1m[22m[30m[47m`inputs`[49m[39m must be a named list!

# use_github_action_pkgcheck [unicode]

    `file_name` must be a character argument

---

    `file_name` must be a single value

---

    `inputs` must be a named list!

# use_github_action_pkgcheck [fancy]

    [1m[22m[30m[47m`file_name`[49m[39m must be a character argument

---

    [1m[22m[30m[47m`file_name`[49m[39m must be a single value

---

    [1m[22m[30m[47m`inputs`[49m[39m must be a named list!

## Checks for [pkgstats (v9.9)](https://github.com/ropensci-review-tools/pkgstats)

git hash: [](https://github.com/ropensci-review-tools/pkgstats/tree/)

- :heavy_check_mark: Package name is available
- :heavy_check_mark: has a 'CITATION' file.
- :heavy_multiplication_x: Package depends on the following obsolete packages: [blah]
- :heavy_multiplication_x: does not have a 'codemeta.json' file.
- :heavy_multiplication_x: does not have a 'contributing' file.
- :heavy_check_mark: uses 'roxygen2'.
- :heavy_check_mark: 'DESCRIPTION' has a URL field.
- :heavy_check_mark: 'DESCRIPTION' has a BugReports field.
- :heavy_multiplication_x: Package has no HTML vignettes
- :heavy_multiplication_x: These functions do not have examples: [pkgstats_from_archive].
- :heavy_check_mark: Package has continuous integration checks.
- :heavy_multiplication_x: Package contains unexpected files.

**Important:** All failing checks above must be addressed prior to proceeding

Package License: GPL-3

---


### 1. Statistical Properties

This package features some noteworthy statistical properties which may need to be clarified by a handling editor prior to progressing.

<details>
<summary>Details of statistical properties (click to open)</summary>
<p>

The package has:

- code in C++ (9% in 3 files) and R (91% in 19 files)
- 1 authors
- no  vignette
- no internal data file
- 9 imported packages
- 11 exported functions (median 43 lines of code)
- 120 non-exported functions in R (median 21 lines of code)
- 12 R functions (median 16 lines of code)

---

Statistical properties of package structure as distributional percentiles in relation to all current CRAN packages
The following terminology is used:
- `loc` = "Lines of Code"
- `fn` = "function"
- `exp`/`not_exp` = exported / not exported

The final measure (`fn_call_network_size`) is the total number of calls between functions (in R), or more abstract relationships between code objects in other languages. Values are flagged as "noteworthy" when they lie in the upper or lower 5th percentile.

|measure                 | value| percentile|noteworthy |
|:-----------------------|-----:|----------:|:----------|
|files_R                 |    19|       79.7|           |
|files_src               |     3|       84.3|           |
|files_vignettes         |     0|        0.0|TRUE       |
|files_tests             |     7|       86.4|           |
|loc_R                   |  2698|       89.0|           |
|loc_src                 |   277|       33.9|           |
|loc_tests               |   266|       61.5|           |
|num_vignettes           |     0|        0.0|TRUE       |
|n_fns_r                 |   131|       82.6|           |
|n_fns_r_exported        |    11|       48.6|           |
|n_fns_r_not_exported    |   120|       87.0|           |
|n_fns_src               |    12|       33.3|           |
|n_fns_per_file_r        |     4|       58.6|           |
|n_fns_per_file_src      |     4|       40.2|           |
|num_params_per_fn       |     1|        1.6|TRUE       |
|loc_per_fn_r            |    23|       66.0|           |
|loc_per_fn_r_exp        |    43|       75.2|           |
|loc_per_fn_r_not_exp    |    22|       66.5|           |
|loc_per_fn_src          |    16|       55.6|           |
|rel_whitespace_R        |    19|       88.9|           |
|rel_whitespace_src      |    24|       41.5|           |
|rel_whitespace_tests    |    27|       64.6|           |
|doclines_per_fn_exp     |    31|       34.8|           |
|doclines_per_fn_not_exp |     0|        0.0|TRUE       |
|fn_call_network_size    |   104|       79.9|           |

---

</p></details>


### 1a. Network visualisation

An interactive visualisation of calls between objects in the package has been uploaded as a workflow artefact. To view it, click on results from the [latest 'pkgcheck' action](network.html), scroll to the bottom, and click on the 'visual-network' artefact.

---

### 2. `goodpractice` and other checks

('goodpractice' not included with these checks)

---

### 3. Other Checks


:heavy_multiplication_x: Package contains the following unexpected files:

- a
- b


:heavy_multiplication_x: Package contains the following (potentially) obsolete packages:

- sp
- rgdal


See our [Recommended Scaffolding](https://devguide.ropensci.org/building.html?q=scaffol#recommended-scaffolding) for alternatives.


---

<details>
<summary>Package Versions</summary>
<p>

|package  |version   |
|:--------|:---------|
|pkgstats |42    |
|pkgcheck |42    |

</p>
</details>
## Checks for [testpkgchecknotapkg (v0.0.0.9000)]()

git hash: [](/tree/)

- :heavy_check_mark: Package name is available
- :heavy_multiplication_x: does not have a 'CITATION' file.
- :heavy_multiplication_x: does not have a 'codemeta.json' file.
- :heavy_multiplication_x: does not have a 'contributing' file.
- :heavy_check_mark: uses 'roxygen2'.
- :heavy_multiplication_x: 'DESCRIPTION' does not have a URL field.
- :heavy_multiplication_x: 'DESCRIPTION' does not have a BugReports field.
- :heavy_multiplication_x: Package has no HTML vignettes
- :heavy_multiplication_x: These functions do not have examples: [test_fn].
- :heavy_multiplication_x: Continuous integration checks unavailable (no URL in 'DESCRIPTION').
- :heavy_multiplication_x: Statistical standards are missing
- :heavy_multiplication_x: This package still has TODO standards and can not be submitted

**Important:** All failing checks above must be addressed prior to proceeding

Package License: GPL-3

---

### 1. rOpenSci Statistical Standards ([`srr` package](https://github.com/ropensci-review-tools/srr))

This package is in the following category:

- *Regression and Supervised Learning*

:heavy_multiplication_x: This package still has TODO standards and can not be submitted

Click to see the [report of author-reported standards compliance of the package with links to associated lines of code](report.html), which can be re-generated locally by running the [`srr_report()` function](https://docs.ropensci.org/srr/reference/srr_report.html) from within a local clone of the repository.

---


### 2. Statistical Properties

This package features some noteworthy statistical properties which may need to be clarified by a handling editor prior to progressing.

<details>
<summary>Details of statistical properties (click to open)</summary>
<p>

The package has:

- code in C++ (72% in 2 files) and R (28% in 4 files)
- 1 authors
- no  vignette
- no internal data file
- 1 imported package
- 1 exported function (median 3 lines of code)
- 2 non-exported functions in R (median 3 lines of code)
- 2 R functions (median 5 lines of code)

---

Statistical properties of package structure as distributional percentiles in relation to all current CRAN packages
The following terminology is used:
- `loc` = "Lines of Code"
- `fn` = "function"
- `exp`/`not_exp` = exported / not exported

The final measure (`fn_call_network_size`) is the total number of calls between functions (in R), or more abstract relationships between code objects in other languages. Values are flagged as "noteworthy" when they lie in the upper or lower 5th percentile.

|measure                 | value| percentile|noteworthy |
|:-----------------------|-----:|----------:|:----------|
|files_R                 |     4|       28.3|           |
|files_src               |     2|       79.1|           |
|files_vignettes         |     0|        0.0|TRUE       |
|files_tests             |     2|       68.6|           |
|loc_R                   |    10|        0.8|TRUE       |
|loc_src                 |    26|        4.0|TRUE       |
|loc_tests               |     6|        4.7|TRUE       |
|num_vignettes           |     0|        0.0|TRUE       |
|n_fns_r                 |     3|        2.5|TRUE       |
|n_fns_r_exported        |     1|        0.0|TRUE       |
|n_fns_r_not_exported    |     2|        2.7|TRUE       |
|n_fns_src               |     2|        4.3|TRUE       |
|n_fns_per_file_r        |     1|        0.2|TRUE       |
|n_fns_per_file_src      |     1|        0.1|TRUE       |
|num_params_per_fn       |     0|        0.0|TRUE       |
|loc_per_fn_r            |     3|        1.1|TRUE       |
|loc_per_fn_r_exp        |     3|        1.5|TRUE       |
|loc_per_fn_r_not_exp    |     3|        1.5|TRUE       |
|loc_per_fn_src          |     5|        5.0|TRUE       |
|rel_whitespace_R        |    40|        4.3|TRUE       |
|rel_whitespace_src      |    27|        5.8|           |
|rel_whitespace_tests    |    17|        2.4|TRUE       |
|doclines_per_fn_exp     |     6|        0.8|TRUE       |
|doclines_per_fn_not_exp |     0|        0.0|TRUE       |
|fn_call_network_size    |     1|       11.4|           |

---

</p></details>


### 2a. Network visualisation

An interactive visualisation of calls between objects in the package has been uploaded as a workflow artefact. To view it, click on results from the [latest 'pkgcheck' action](network.html), scroll to the bottom, and click on the 'visual-network' artefact.

---

### 3. `goodpractice` and other checks

('goodpractice' not included with these checks)

---

<details>
<summary>Package Versions</summary>
<p>

|package  |version   |
|:--------|:---------|
|pkgstats |42    |
|pkgcheck |42    |
|srr      |42    |

</p>
</details>
---
title: "pgkcheck"
output:
  md_document:
    variant: gfm

  rmarkdown::html_vignette:
    self_contained: no
---

# pkgcheck

<!-- badges: start -->
[![R build status](https://github.com/ropensci-review-tools/pkgcheck/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci-review-tools/pkgcheck/actions?query=workflow%3AR-CMD-check)
[![gitlab push](https://github.com/ropensci-review-tools/pkgcheck/workflows/push-to-gitlab/badge.svg)](https://github.com/ropensci-review-tools/pkgcheck/actions?query=workflow%3Apush-to-gitlab)
[![codecov](https://codecov.io/gh/ropensci-review-tools/pkgcheck/branch/main/graph/badge.svg)](https://codecov.io/gh/ropensci-review-tools/pkgcheck)
[![Project Status: Concept](https://www.repostatus.org/badges/latest/concept.svg)](https://www.repostatus.org/#concept)
<!-- badges: end -->

Check whether a package is ready for submission to
[rOpenSci](https://ropensci.org)'s peer review system. The primary function
collates the output of
[`goodpractice`](https://github.com/mangothecat/goodpractice), including `R CMD
check` results, a number of statistics via the [`pkgstats`
package](https://github.com/ropensci-review-tools/pkgstats), and checks for
package structure expected for rOpenSci submissions. The output of this
function immediately indicates whether or not a package is "Ready to Submit".


## Installation

The easiest way to install this package is via the [associated
`r-universe`](https://ropensci-review-tools.r-universe.dev/ui#builds). As
shown there, simply enable the universe with

```{r options, eval = FALSE}
options (repos = c (
    ropenscireviewtools = "https://ropensci-review-tools.r-universe.dev",
    CRAN = "https://cloud.r-project.org"
))
```

And then install the usual way with,

```{r install, eval = FALSE}
install.packages ("pkgcheck")
```

Alternatively, the package can be installed by running one of the following
lines:



```{r remotes, eval = FALSE}
remotes::install_github ("ropensci-review-tools/pkgcheck")
pak::pkg_install ("ropensci-review-tools/pkgcheck")
```

The package can then loaded for use with

```{r library, eval = FALSE}
library (pkgcheck)
```


## Setup

The [`pkgstats` package](https://github.com/ropensci-review-tools/pkgstats)
also requires the system libraries [`ctags`](https://ctags.io) and [GNU
`global`](https://www.gnu.org/software/global/) to be installed. Procedures to
install these libraries on various operating systems are described in the
[`pkgstats` README](https://docs.ropensci.org/pkgstats). This package also uses
the [GitHub GraphQL API](https://developer.github.com/v4) which requires
a local GitHub token to be stored with an unambiguous name including `GITHUB`,
such as `GITHUB_TOKEN` (recommended), or `GITHUB_PAT` (for Personal
Authorization Token). This can be obtained from GitHub (via your user
settings), and stored using

```{r ghtok, eval = FALSE}
Sys.setenv ("GITHUB_TOKEN" = "<my_token>")
```

This can also be set permanently by putting this line in your `~/.Renviron`
file (or creating this if it does not yet exist). Once `pkgstats` has been
successfully installed, the `pkgcheck` package can then be loaded via
a `library` call:

```{r load, echo = FALSE, message = FALSE}
devtools::load_all (".", export_all = FALSE)
```
```{r load-fakey, eval = FALSE}
library (pkgcheck)
```


## Usage

The package primarily has one function, `pkgcheck`, which accepts the single
argument, `path`, specifying the local location of a git repository to be
analysed. The following code generates a reproducible report by first
downloading a local clone of a repository called
[`srr-demo`](https://github.com/mpadge/srr-demo), which contains the skeleton
of an [`srr` (Software Review Roclets)
package](https://github.com/ropensci-review-tools/srr), generated with the
[`srr_stats_pkg_skeleton()`
function](https://docs.ropensci.org/srr/reference/srr_stats_pkg_skeleton.html):

```{r pkgcheck, echo = TRUE, message = FALSE, include = FALSE}
mydir <- file.path (tempdir (), "srr-demo")
gert::git_clone ("https://github.com/mpadge/srr-demo", path = mydir)
x <- pkgcheck (mydir)
```

That object has default `print` and `summary` methods. The latter can be used
to simply check whether a package is ready for submission:

```{r summary, collapse = TRUE}
summary (x)
```

A package may only be submitted when the summary contains all ticks and no
cross symbols. (These symbols are colour-coded with green ticks and red crosses
when generated in a terminal; GitHub markdown only renders them in
black-and-white.) The object returned from the `pkgcheck` function is a complex
nested list with around a dozen primary components. Full information can be
obtained by simply calling the default `print` method by typing the object name
(`x`).

## What is checked?

### Summary of Check Results

Calling `summary()` on the object returned by the [`pkgcheck()`
function](https://docs.ropensci.org/pkgcheck/reference/pkgcheck.html) will
generate a checklist like that shown above. This checklist will also be
automatically generated when a package is first submitted to rOpenSci, and is
used by the editors to assess whether to process a submission. Authors must
ensure prior to submission that there are no red crosses in the resultant list.
(In the unlikely circumstances that a package is unable to pass particular
checks, explanations should be given upon submission about why those checks
fail, and why review may proceed in spite of such failures.)

The full list of checks which packages are expected to pass currently includes:

1. Package must use [`roxygen2`](https://devguide.ropensci.org/building.html#roxygen2-use) for documentation.
2. Package must have a [`contributing.md`
   file](https://devguide.ropensci.org/collaboration.html#contributing-guide).
3. Package must have a [`CITATION` file in the `inst`
   directory](https://cran.r-project.org/doc/manuals/r-release/R-exts.html#CITATION-files).
4. Package must have a [`codemeta.json`
   file](https://devguide.ropensci.org/building.html#creating-metadata-for-your-package).
5. All exported functions must include examples in their documentation.
6. Left-assign operators must be used consistently throughout all code (so
   either all `=` or all `<-`, but not a mixture of both).
7. Package `DESCRIPTION` file must have a "URL" field.
8. Package `DESCRIPTION` file must have a "BugReports" field.
9. Package name must be available (or package must already be) on CRAN.
10. Package must have continuous integration tests.
11. Package must have test coverage of at least 75%.
12. No function should have a [cyclomatic
    complexity](https://github.com/MangoTheCat/cyclocomp) of 15 or greater.
13. `R CMD check` (implemented via the [`rcmdcheck`
    package](https://r-lib.github.io/rcmdcheck/)) must generate no warnings or
    errors.
14. All statistical standards must be documented, as confirmed by the
    [`srr::srr_pre_submit()`
    function](https://docs.ropensci.org/srr/reference/srr_stats_pre_submit.html).

Several of these are by default only shown when they fail; absence from
a resultant checklist may be interpreted to indicate successful checks.

### Detailed Check Results

Full details of check results can be seen by `print`-ing the object returned by
the [`pkgcheck()`
function](https://docs.ropensci.org/pkgcheck/reference/pkgcheck.html) (or just
by typing the name of this object in the console.) This object is itself a list
including the following items:

```{r check-items}
names (x)
```

The first four of these contain information on the package. The remainder include:

- `summary` containing summary statistics on package structure (such as lines
  of code, and numbers of files and functions).
- `git` containing information and statistics on the `git` repository
  associated with a package (if any), including date of first commit, total
  number of commits and of committers.
- `srr` summarising information on statistical packages generated by the [`srr`
  package](https://github.com/ropensci-review-tools/srr), including a link to
  a locally-generated HTML report on standards compliance.
- `file_list` as a list of binary flags denoting presence of absence of all
  files required for an R package to be submitted to rOpenSci
- `fns_have_exs` with a list of any functions which do not have documented
  examples
- `left_assigns` with a logical value indicating whether or not global
  assignment operators (`<<-`) are used; and a tally of the two types of
  left-assignment operators (`<-` and `=`). One one of these should be
  non-zero, reflecting consistent usage of the same type.
- `pkgstats` containing several statistics of package structure, along with
  associated percentiles assessed against the entire distribution of all
  current CRAN packages.
- `network_file` with the path to a local HTML file containing
  a [.visjs](https://visjs.org/) representation of the network of relationships
  between objects (such as functions) within a package, within and between each
  computer language used in the package.
- `gp` containing the output of the [`goodpractice`
  pacakge](http://mangothecat.github.io/goodpractice/), itself including
  results from:
    - the [`rcmdcheck` pacakge](https://r-lib.github.io/rcmdcheck) for running
      `R CMD check`
    - the [`covr` package](https://covr.r-lib.org) for assessing code coverage
    - the [`cyclocomp` package](https://github.com/MangoTheCat/cyclocomp) for
      quantifying the cyclomatic complexity of package functions
    - the [`desc` package](https://github.com/r-lib/desc#readme) for analysing
      the structure of `DESCRIPTION` files
    - the [`lintr` package](https://github.com/jimhester/lintr)
    - additional checks applied by the [`goodpractice`
      pacakge](http://mangothecat.github.io/goodpractice/)

Note that results from [`lintr` package](https://github.com/jimhester/lintr)
are **not** reported in the check summary, and that [`lintr`
results](https://github.com/jimhester/lintr) are reported only in the detailed
results, and have no influence on whether a package passes the summary checks.


### HTML-formatted check results

Printing `pkgcheck` results to screen is nice, but sometimes it's even better
to have a nicely formatted full-screen representation of check results. The
package includes an additional function,
[`checks_to_markdown()`](https://docs.ropensci.org/pkgcheck/reference/checks_to_markdown.html),
with a parameter, `render`, which can be set to `TRUE` to automatically render
a HTML-formatted representation of the check results, and open it in a
browser. The formatting differs only slightly from the terminal output, mostly
through the components of
[`goodpractice`](http://mangothecat.github.io/goodpractice/) being divided into
distinct headings, with explicit statements in cases where components pass all
checks (the default screen output inherits directly from that package, and only
reports components which *do not* pass all checks).

This
[`checks_to_markdown()`](https://docs.ropensci.org/pkgcheck/reference/checks_to_markdown.html)
function returns the report in markdown format, suitable for pasting directly
into a GitHub issue, or other markdown-compatible place. (The [`clipr`
package](https://github.com/mdlincoln/clipr) can be used to copy this directly
to your local clipboard with `write_clip(md)`, where `md` is the output of
`checks_to_markdown()`.)


## Caching and running `pkgcheck` in the background

Running the [`pgkcheck`
function](https://docs.ropensci.org/pkgcheck/reference/pkgcheck.html) can
be time-consuming, primarily because the
[`goodpractice`](https://github.com/mangothecat/goodpractice) component runs
both a full `R CMD check`, and calculates code coverage of all tests. To avoid
re-generating these results each time, the package saves previous reports to
a local cache, in a `pkgcheck` subdirectory of the location determined by

```{r, rappdirs, eval = FALSE}
rappdirs::user_cache_dir ()
```

As explained in the help file for that function, these locations are:

System | location
--- | ---
Mac OS X | `~/Library/Caches/pkgcheck`
Linux | `~/.cache/pkgcheck`
Win XP | `C:\\Documents and Settings\\<username>\\Local Settings\\Application Data\\<AppAuthor>\\pkgcheck\\Cache`
Vista | `C:\\Users\\<username>\\AppData\\Local\\<AppAuthor>\\pkgcheck\\Cache`


You may manually erase the contents of this `pkgcheck` subdirectory at any time
at no risk beyond additional time required to re-generate contents. This
default location may also be over-ridden by setting an environmental variable
named `PKGCHECK_CACHE_DIR`. By default checks presume packages use `git` for
version control, with checks updated only when code is updated via `git
commit`. Checks for packages that do not use `git` are updated when any files
are modified.

The first time
[`pkgcheck()`](https://docs.ropensci.org/pkgcheck/reference/pkgcheck.html)
is applied to a package, the checks will be stored in the cache directory.
Calling that function a second time will then load the cached results, and so
enable checks to be returned much faster. For code which is frequently updated,
such as for packages working on the final stages prior to submission, it may
still be necessary to repeatedly call
[`pkgcheck()`](https://docs.ropensci.org/pkgcheck/reference/pkgcheck.html)
after each modification, a step which may still be inconveniently
time-consuming. To facilitate frequent re-checking, the package also has
a [`pkgcheck_bg()`
function](https://docs.ropensci.org/pkgcheck/reference/pkgcheck_bg.html)
which is effectively identical to the main 
[`pkgcheck()`
function](https://docs.ropensci.org/pkgcheck/reference/pkgcheck.html),
except it runs in the background, enabling you to continue coding while checks
are running.

The [`pkgcheck_bg()`
function](https://docs.ropensci.org/pkgcheck/reference/pkgcheck_bg.html)
returns a handle to the [`callr::r_bg()`
process](https://callr.r-lib.org/reference/r_bg.html) in which the checks are
running. Typing the name of the returned object will immediately indicate
whether the checks are still running, or whether they have finished. That
handle is itself an [`R6` object](http://r6.r-lib.org/) with a number of
methods, notably including
[`get_result()`](https://callr.r-lib.org/reference/get_result.html) which can
be used to access the checks once the process has finished. Alternatively, as
soon as the background process, the normal (foreground) 
[`pkgcheck()`
function](https://docs.ropensci.org/pkgcheck/reference/pkgcheck.html) may
be called to quickly re-load the cached results.

## Prior Work

[The `checklist` package](https://github.com/inbo/checklist) for "checking
packages and R code".

## Code of Conduct

Please note that this package is released with a [Contributor Code of
Conduct](https://ropensci.org/code-of-conduct/). By contributing to this
project, you agree to abide by its terms.

## Contributors


<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->

All contributions to this project are gratefully acknowledged using the [`allcontributors` package](https://github.com/ropenscilabs/allcontributors) following the [all-contributors](https://allcontributors.org) specification. Contributions of any kind are welcome!

### Code

<table>

<tr>
<td align="center">
<a href="https://github.com/mpadge">
<img src="https://avatars.githubusercontent.com/u/6697851?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci-review-tools/pkgcheck/commits?author=mpadge">mpadge</a>
</td>
<td align="center">
<a href="https://github.com/maelle">
<img src="https://avatars.githubusercontent.com/u/8360597?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci-review-tools/pkgcheck/commits?author=maelle">maelle</a>
</td>
<td align="center">
<a href="https://github.com/assignUser">
<img src="https://avatars.githubusercontent.com/u/16141871?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci-review-tools/pkgcheck/commits?author=assignUser">assignUser</a>
</td>
<td align="center">
<a href="https://github.com/annakrystalli">
<img src="https://avatars.githubusercontent.com/u/5583057?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci-review-tools/pkgcheck/commits?author=annakrystalli">annakrystalli</a>
</td>
<td align="center">
<a href="https://github.com/noamross">
<img src="https://avatars.githubusercontent.com/u/571752?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci-review-tools/pkgcheck/commits?author=noamross">noamross</a>
</td>
</tr>

</table>


### Issue Authors

<table>

<tr>
<td align="center">
<a href="https://github.com/piyalkarum">
<img src="https://avatars.githubusercontent.com/u/48254643?u=370433a2ace6a030f2551575bc08fa53664fbd8f&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci-review-tools/pkgcheck/issues?q=is%3Aissue+author%3Apiyalkarum">piyalkarum</a>
</td>
<td align="center">
<a href="https://github.com/christophsax">
<img src="https://avatars.githubusercontent.com/u/1390827?u=ce6363f6da758d1bb85987d021cacc34a81c8837&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci-review-tools/pkgcheck/issues?q=is%3Aissue+author%3Achristophsax">christophsax</a>
</td>
<td align="center">
<a href="https://github.com/steffilazerte">
<img src="https://avatars.githubusercontent.com/u/14676081?u=0c1a4469750eba5a6f398c35ba34a2dd600a873e&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci-review-tools/pkgcheck/issues?q=is%3Aissue+author%3Asteffilazerte">steffilazerte</a>
</td>
</tr>

</table>


### Issue Contributors

<table>

<tr>
<td align="center">
<a href="https://github.com/dgkf">
<img src="https://avatars.githubusercontent.com/u/18220321?u=bef717254e5b877159fa712e2b8ad6952c816064&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci-review-tools/pkgcheck/issues?q=is%3Aissue+commenter%3Adgkf">dgkf</a>
</td>
<td align="center">
<a href="https://github.com/cboettig">
<img src="https://avatars.githubusercontent.com/u/222586?u=dfbe54d3b4d538dc2a8c276bb5545fdf4684752f&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci-review-tools/pkgcheck/issues?q=is%3Aissue+commenter%3Acboettig">cboettig</a>
</td>
<td align="center">
<a href="https://github.com/jhollist">
<img src="https://avatars.githubusercontent.com/u/5438539?u=815aa29d708acd16180967c0ffaf81fc64a08bf4&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci-review-tools/pkgcheck/issues?q=is%3Aissue+commenter%3Ajhollist">jhollist</a>
</td>
</tr>

</table>

<!-- markdownlint-enable -->
<!-- prettier-ignore-end -->
<!-- ALL-CONTRIBUTORS-LIST:END -->
---
title: "Extending or modifying checks"
author: 
  - "Mark Padgham"
date: "`r Sys.Date()`"
vignette: >
  %\VignetteIndexEntry{How to extend checks}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set (
    collapse = TRUE,
    warning = TRUE,
    message = TRUE,
    width = 120,
    comment = "#>",
    fig.retina = 2,
    fig.path = "README-"
)
options (repos = c (
    ropenscireviewtools = "https://ropensci-review-tools.r-universe.dev",
    CRAN = "https://cloud.r-project.org"
))
```

This vignette describes how to modify or extend the existing suite of checks
implemented by `pkgcheck`. Each of the internal checks is defined in a separate
file in the `R` directory of this package with the prefix of `check_` (or
`checks_` for files which define multiple, related checks). Checks only require
two main functions, the first defining the check itself, and the second
defining `summary` and `print` methods based on the result of the first
function. The check functions must have a prefix `pkgchk_`, and the second
functions defining output methods specifying must have a prefix
`output_pkgchk_`. These two kind of function are now described in the following
two sections. 

Both of these functions must also accept a single input parameter of a
`pkgcheck` object, by convention named `checks`. This object is a list of four
main items:

1. `pkg` which summarises data extracted from
   [`pkgstats::pkgstats()`](https://docs.ropensci.org/pkgstats/reference/pkgstats.html),
   and includes essential information on the package being checked.
2. `info` which contains information used in checks, including `info$git`
   detailing git repository information, `info$pkgstats` containing a summary
   of a few statistics generated from
   [`pkgstats::pkgstats()`](https://docs.ropensci.org/pkgstats/reference/pkgstats.html),
   along with statistical comparisons against distributions from all current
   CRAN packages, an `info$network_file` specifying a local directory to a
   [`vis.js`](https://visjs.org) visualisation of the function call network of
   the package, and an `info$badges` item containing information from GitHub
   workflows and associated badges, where available.
3. `checks` which contains a list of all objects returned from all
   `pkgchk_...()` functions, which are used as input to `output_pkgchk_...()`
   functions.
4. `meta` containing a named character vector of versions of the core packages
   used in `pkgcheck`.

`pkgcheck` objects generally also include a fifth item, `goodpractice`,
containing the results of [`goodpractice`
checks](https://github.com/MangoTheCat/goodpractice). The `checks` item passed
to each `pkgchk_...()` function contains all information on the `package`,
`info`, `meta`, and (optionally) `goodpractice` items. Checks may use any of
this information, or even add additional information as demonstrated below. The
`checks$checks` list represents the output of check functions, and may not be
used in any way within `pkgchk_...()` functions.

<details>
<summary>Click here to see structure of full `pkgcheck` object</summary>
<p>


This is the output of applying `pkgcheck` to a package generated with the
[`srr` function
`srr_stats_pkg_skeleton()`](https://docs.ropensci.org/srr/reference/srr_stats_pkg_skeleton.html),
with `goodpractice = FALSE` to suppress that part of the results.

```{r check_file, echo = FALSE}
here <- rprojroot::find_root (rprojroot::is_r_package)
check_file <- file.path (here, "vignettes", "checks.Rds")
```


```{r pkgcheck-dummy-data, echo = FALSE, eval = !file.exists (check_file)}
d <- srr::srr_stats_pkg_skeleton (pkg_name = "dummypkg")
roxygen2::roxygenise (d)
checks <- pkgcheck::pkgcheck (d, goodpractice = FALSE)
saveRDS (checks, check_file)
```
                       
```{r pkgcheck-str, echo = FALSE, cache = FALSE}
print (str (readRDS (check_file)))
```

</p></details>




## 1. The check function



An example is the
check for whether a package has a citation, [defined in
`R/check_has_citation.R`](https://github.com/ropensci-review-tools/pkgcheck/blob/main/R/check-has-citation.R):

```{r, cache = FALSE, echo = FALSE}
knitr::read_chunk ("../R/check-has-citation.R")
knitr::read_chunk ("../R/check-scrap.R")
```
```{r pkgchk-citation}
```

This check is particularly simple, because a `"CITATION"` file [must have
exactly that name, and must be in the `inst`
sub-directory](https://cran.r-project.org/doc/manuals/R-exts.html#CITATION-files).
This function returns a simple logical of `TRUE` if the expected `"CITATION"`
file is present, otherwise it returns `FALSE`. This function, and all functions
beginning with the prefix `pkgchk_`, will be automatically called by the main
`pkgcheck()` function, and the value stored in `checks$checks$has_citation`.
The name of the item within the `checks$checks` list is the name of the
function with the `pkgchk_` prefix removed.

A more complicated example is the function to check whether a package contains
files which should not be there -- internally called "scrap" files. The check
function itself, [defined in
`R/check-scrap.R`](https://github.com/ropensci-review-tools/pkgcheck/blob/main/R/check-scrap.R),
checks for the presence of files matching an internally-defined list including
files used to locally cache folder thumbnails such as `".DS_Store"` or
`"Thumbs.db"`. The function returns a character vector of the names of any
"scrap" files which can be used by the `print` method to provide details of
files which should be removed. This illustrates the first general principle of
these check functions; that,

::: {.alert .alert-info}
- *Any information needed when summarising or printing the check result should
  be returned from the main check function.*
:::


A second important principle is that,

::: {.alert .alert-info}
- *Check functions should never return `NULL`, rather should always return an
  empty vector (such as `integer(0)`)*.
:::

The following section considers how these return values from check functions
are converted to `summary` and `print` output.

## 2. The output function

All `output_pkgchk_...()` functions must also accept the single input parameter
of `checks`, in which the `checks$checks` sub-list will already have been
populated by calling all `pkgchk_...()` functions described in the previous
section. The `pkgchk_has_citation()` function will create an entry of
`checks$checks$has_citation` which contains the binary flag indicating whether
or not a `"CITATION"` file is present. Similarly, the [the `pkgchk_has_scrap()`
function](https://github.com/ropensci-review-tools/pkgcheck/blob/main/R/check-scrap.R)
will create `checks$checks$has_scrap` which will contain names of any scrap
files present, and a length-zero vector otherwise.

The `output_pkgchk_has_citation()` function then looks like this:

```{r output-pkgchk-citation}
```

The first lines are common to all `output_pkgchk_...()` functions, and define
the generic return object. This object must be a list with the following three
items:

1. `check_pass` as binary flag indicating whether or not a check was passed;
2. `summary` containing text used to generate the `summary` output; and
3. `print` containing information used to generate the `print` output, itself a
   `list` of the following items:
    - A `msg_pre` to display at the start of the `print` result;
    - An `object` to be printed, such as a vector of values, or a `data.frame`.
    - A `msg_post` to display at the end of the `print` result following the
      `object`.

`summary` and `print` methods may be suppressed by assigning values of `""`.
The above example of `pkgcheck_has_citation` has `print = ""`, and so no
information from this check will appear as output of the `print` method. The
`summary` field has a value that is specified for both `TRUE` and `FALSE`
values of `check_pass`. The value is determined by the result of the main
`pkgchk_has_citation()` call, and is converted into a green tick if `TRUE`, or
a red cross if `FALSE`.

Checks for which `print` information is desired require a non-empty `print`
item, as in the [`output_pkgchk_has_scrap()`
function](https://github.com/ropensci-review-tools/pkgcheck/blob/main/R/check-scrap.R):

```{r output-pkgchk-scrap}
```

In this case, both `summary` and `print` methods are only triggered `if
(!out$check_pass)` -- so only if the check fails. The `print` method generates
the heading specified in `out$print$msg_pre`, with any vector-valued objects
stored in the corresponding `obj` list item displayed as formatted lists.
A package with "scrap" files, `"a"` and `"b"`, would thus have `out$print$obj
<- c ("a", "b")`, and when printed would look like this:

```{r scrap-out, echo = FALSE}
cli::cli_alert_danger ("Package contains the following unexpected files:")
cli::cli_ul ()
cli::cli_li (c ("a", "b"))
cli::cli_end ()
```

This formatting is also translated into corresponding markdown and HTML
formatting in [the `checks_to_markdown()`
function](https://github.com/ropensci-review-tools/pkgcheck/blob/main/R/format-checks.R).

The design of these `pkgchk_` and `output_pkgchk_` functions aims to make the
package readily extensible, and we welcome discussions about developing new
checks. The primary criterion for new package-internal checks is that they must
be of very general applicability, in that they should check for a condition
that *almost* every package should or should not meet.

The package also has a mechanism to easily incorporate more specific,
locally-defined checks, as explored in the following section.

## 3. Creating new checks

(Coming soon ...)
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{list_pkgchecks}
\alias{list_pkgchecks}
\title{List all checks currently implemented}
\usage{
list_pkgchecks(quiet = FALSE)
}
\arguments{
\item{quiet}{If \code{TRUE}, print all checks to screen. Function invisibly
returns list of checks regardless.}
}
\value{
Character vector of names of all checks (invisibly)
}
\description{
List all checks currently implemented
}
\seealso{
Other extra: 
\code{\link{checks_to_markdown}()},
\code{\link{logfile_names}()},
\code{\link{render_markdown}()}
}
\concept{extra}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/github.R
\name{get_latest_commit}
\alias{get_latest_commit}
\title{get_latest_commit}
\usage{
get_latest_commit(org, repo)
}
\arguments{
\item{org}{Github organization}

\item{repo}{Github repository}
}
\value{
Details of latest commit including OID hash
}
\description{
get_latest_commit
}
\note{
This returns the latest commit from the default branch as specified on
GitHub, which will not necessarily be the same as information returned from
\code{gert::git_info} if the \code{HEAD} of a local repository does not point to the
same default branch.
}
\seealso{
Other github: 
\code{\link{get_default_branch}()},
\code{\link{get_gh_token}()},
\code{\link{use_github_action_pkgcheck}()}
}
\concept{github}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format-checks.R
\name{checks_to_markdown}
\alias{checks_to_markdown}
\title{Convert checks to markdown-formatted report}
\usage{
checks_to_markdown(checks, render = FALSE)
}
\arguments{
\item{checks}{Result of main \link{pkgcheck} function}

\item{render}{If \code{TRUE}, render output as \code{html} document and open in
browser.}
}
\value{
Markdown-formatted version of check report
}
\description{
Convert checks to markdown-formatted report
}
\seealso{
Other extra: 
\code{\link{list_pkgchecks}()},
\code{\link{logfile_names}()},
\code{\link{render_markdown}()}
}
\concept{extra}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/github.R
\name{get_gh_token}
\alias{get_gh_token}
\title{Get GitHub token}
\usage{
get_gh_token(token_name = "")
}
\arguments{
\item{token_name}{Optional name of token to use}
}
\value{
The value of the GitHub access token extracted from environment
variables.
}
\description{
Get GitHub token
}
\seealso{
Other github: 
\code{\link{get_default_branch}()},
\code{\link{get_latest_commit}()},
\code{\link{use_github_action_pkgcheck}()}
}
\concept{github}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pkgcheck-bg.R
\name{pkgcheck_bg}
\alias{pkgcheck_bg}
\title{Generate report on package compliance with rOpenSci Statistical Software
requirements as background process}
\usage{
pkgcheck_bg(path)
}
\arguments{
\item{path}{Path to local repository}
}
\value{
A \pkg{processx} object connecting to the background process
generating the main \link{pkgcheck} results (see Note).
}
\description{
Generate report on package compliance with rOpenSci Statistical Software
requirements as background process
}
\note{
The return object will by default display whether it is still running,
or whether it has finished. Once it has finished, the results can be obtained
by calling \verb{$get_result()}, or the main \link{pkgcheck} function can be
called to quickly retrieve the main results from local cache.

This function does not accept the \code{extra_env} parameter of the main
\link{pkgcheck} function, and can not be used to run extra, locally-defined
checks.
}
\seealso{
Other pkgcheck_fns: 
\code{\link{pkgcheck}()}
}
\concept{pkgcheck_fns}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/github.R
\name{use_github_action_pkgcheck}
\alias{use_github_action_pkgcheck}
\title{Use pkgcheck Github Action}
\usage{
use_github_action_pkgcheck(
  dir = ".github/workflows",
  overwrite = FALSE,
  file_name = "pkgcheck.yaml",
  inputs = NULL
)
}
\arguments{
\item{dir}{Directory the file is written to.}

\item{overwrite}{Overwrite existing file?}

\item{file_name}{Name of the workflow file.}

\item{inputs}{Named list of inputs to the
\code{ropensci-review-tools/pkgcheck-action}. See details below.}
}
\value{
The path to the new file, invisibly.
}
\description{
Creates a Github workflow file in \code{dir} integrate \code{\link[=pkgcheck]{pkgcheck()}} into your CI.
}
\details{
For more information on the action and advanced usage visit the
action
\href{https://github.com/ropensci-review-tools/pkgcheck-action}{repository}.
}
\section{Inputs}{

Inputs with description and default values. Pass all values as strings, see
examples.\if{html}{\out{<div class="sourceCode yaml">}}\preformatted{inputs:
  ref:
    description: "The ref to checkout and check. Set to empty string to skip
     checkout."
    default: "$\{\{ github.ref \}\}"
    required: true
  post-to-issue:
    description: "Should the pkgcheck results be posted as an issue?"
    # If you use the 'pull_request' trigger and the PR is from outside the
    # repo
    # (e.g. a fork), the job will fail due to permission issues
    # if this is set to 'true'. The default will prevent this.
    default: $\{\{ github.event_name != 'pull_request' \}\}
    required: true
  issue-title:
    description: "Name for the issue containing the pkgcheck results. Will
    be created or updated."
    # This will create a new issue for every branch, set it to something
    # fixed
    # to only create one issue.
    default: "pkgcheck results - $\{\{ github.ref_name \}\}"
    required: true
  summary-only:
    description: "Only post the check summary to issue. Set to false to get
     the full results in the issue."
    default: true
    required: true
}\if{html}{\out{</div>}}
}

\examples{
\dontrun{
use_github_action_pkgcheck (inputs = list (`post-to-issue` = "false"))
}
}
\seealso{
Other github: 
\code{\link{get_default_branch}()},
\code{\link{get_gh_token}()},
\code{\link{get_latest_commit}()}
}
\concept{github}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cache.R
\name{logfile_names}
\alias{logfile_names}
\title{Set up stdout & stderr cache files for \code{r_bg} process}
\usage{
logfile_names(path)
}
\arguments{
\item{path}{Path to local repository}
}
\value{
Vector of two strings holding respective local paths to \code{stdout} and
\code{stderr} files for \code{r_bg} process controlling the main \link{pkgcheck}
function when executed in background mode.
}
\description{
Set up stdout & stderr cache files for \code{r_bg} process
}
\note{
These files are needed for the \pkg{callr} \code{r_bg} process which
controls the main \link{pkgcheck}. The \code{stdout} and \code{stderr} pipes from the
process are stored in the cache directory so they can be inspected via their
own distinct endpoint calls.
}
\seealso{
Other extra: 
\code{\link{checks_to_markdown}()},
\code{\link{list_pkgchecks}()},
\code{\link{render_markdown}()}
}
\concept{extra}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pkgcheck-fn.R
\name{pkgcheck}
\alias{pkgcheck}
\title{Generate report on package compliance with rOpenSci Statistical Software
requirements}
\usage{
pkgcheck(path = ".", goodpractice = TRUE, extra_env = .GlobalEnv)
}
\arguments{
\item{path}{Path to local repository}

\item{goodpractice}{If \code{FALSE}, skip goodpractice checks. May be useful in
development stages to more quickly check other aspects.}

\item{extra_env}{Additional environments from which to collate checks. Other
package names may be appended using \code{c}, as in \code{c(.GlobalEnv, "mypkg")}.}
}
\value{
A \code{pkgcheck} object detailing all package assessments automatically
applied to packages submitted for peer review.
}
\description{
Generate report on package compliance with rOpenSci Statistical Software
requirements
}
\seealso{
Other pkgcheck_fns: 
\code{\link{pkgcheck_bg}()}
}
\concept{pkgcheck_fns}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/github.R
\name{get_default_branch}
\alias{get_default_branch}
\title{get_default_branch}
\usage{
get_default_branch(org, repo)
}
\arguments{
\item{org}{Github organization}

\item{repo}{Github repository}
}
\value{
Name of default branch on GitHub
}
\description{
get_default_branch
}
\note{
This function is not intended to be called directly, and is only
exported to enable it to be used within the \pkg{plumber} API.
}
\seealso{
Other github: 
\code{\link{get_gh_token}()},
\code{\link{get_latest_commit}()},
\code{\link{use_github_action_pkgcheck}()}
}
\concept{github}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format-checks.R
\name{render_markdown}
\alias{render_markdown}
\title{render markdown-formatted input into 'html'}
\usage{
render_markdown(md, open = TRUE)
}
\arguments{
\item{md}{Result of \link{checks_to_markdown} function.}

\item{open}{If \code{TRUE}, open \code{hmtl}-rendered version in web browser.}
}
\value{
(invisible) Location of \code{.html}-formatted version of input.
}
\description{
render markdown-formatted input into 'html'
}
\seealso{
Other extra: 
\code{\link{checks_to_markdown}()},
\code{\link{list_pkgchecks}()},
\code{\link{logfile_names}()}
}
\concept{extra}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pkgcheck-package.R
\docType{data}
\name{pkgstats_data}
\alias{pkgstats_data}
\title{Statistics on all CRAN packages from 'pkgstats'}
\format{
An object of class \code{grouped_df} (inherits from \code{tbl_df}, \code{tbl}, \code{data.frame}) with 21091 rows and 91 columns.
}
\usage{
pkgstats_data
}
\description{
Statistics on all CRAN packages from 'pkgstats'
}
\concept{data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pkgcheck-package.R
\docType{package}
\name{pkgcheck-package}
\alias{pkgcheck-package}
\alias{_PACKAGE}
\title{pkgcheck: rOpenSci Package Checks}
\description{
Check whether a package is ready for submission to rOpenSci’s peer review system.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://docs.ropensci.org/pkgcheck/}
  \item \url{https://github.com/ropensci-review-tools/pkgcheck}
  \item Report bugs at \url{https://github.com/ropensci-review-tools/pkgcheck/issues}
}

}
\author{
\strong{Maintainer}: Mark Padgham \email{mark.padgham@email.com} (\href{https://orcid.org/0000-0003-2172-5265}{ORCID})

Authors:
\itemize{
  \item Maëlle Salmon
  \item Jacob Wujciak-Jens \email{jacob@wujciak.de} (\href{https://orcid.org/0000-0002-7281-3989}{ORCID})
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read-books.R
\name{read_pkg_guide}
\alias{read_pkg_guide}
\title{Browse packaging guidelines}
\usage{
read_pkg_guide(which = c("release", "dev"))
}
\arguments{
\item{which}{Whether to read the released or "dev" development version.}
}
\description{
Browse packaging guidelines
}
