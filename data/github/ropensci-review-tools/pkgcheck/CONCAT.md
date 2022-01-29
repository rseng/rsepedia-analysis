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
