
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `rromeo` – an R interface for SHERPA/RoMEO API

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Travis build
status](https://travis-ci.org/ropensci/rromeo.svg?branch=master)](https://travis-ci.org/ropensci/rromeo)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/ropensci/rromeo?branch=master&svg=true)](https://ci.appveyor.com/project/ropensci/rromeo)
[![codecov](https://codecov.io/gh/ropensci/rromeo/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/rromeo)
[![cran
checks](https://cranchecks.info/badges/summary/rromeo)](https://cran.r-project.org/web/checks/check_results_rromeo.html)
[![CRAN-version](https://www.r-pkg.org/badges/version/rromeo)](https://cran.r-project.org/package=rromeo)
[![](https://badges.ropensci.org/285_status.svg)](https://github.com/ropensci/onboarding/issues/285)

`rromeo` is an R client for the [SHERPA/RoMEO
API](http://www.sherpa.ac.uk/romeo/index.php?la=en&fIDnum=&mode=simple).
SHERPA/RoMEO is a database that gives information on editorial policies
of scientific journals regarding the archival of preprint, postprint and
publishers’ manuscripts. `rromeo` is aimed at scientists interested in
archival practices of scientific journals, such as professionals of
[scientometrics](https://en.wikipedia.org/wiki/Scientometrics) but also
at scientist of specific fields interested in the practices of their
fields.

## Install

The latest stable release of `rromeo` is available on CRAN and can be
installed with:

``` r
install.packages("rromeo")
```

You can also install the development version from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("ropensci/rromeo")
```

## API Key

Note that SHERPA/RoMEO lets you run 500 requests per day per IP address,
by [registering for a free API
key](http://www.sherpa.ac.uk/romeo/apiregistry.php) you can bypass this
limit.

`rromeo` can use your registered SHERPA/RoMEO API key; you can either
pass it as a string when querying the data with the argument `key`:

``` r
rr_journal_name("Journal of Geology", key = "Iq83AIL5bss")
```

or you can specify the environment variable `SHERPAROMEO_KEY` in an
`.Rprofile` or in an `.Renviron` file and `rromeo` will automatically
retrieve the API key. See the [specific
vignette](https://docs.ropensci.org/rromeo/articles/setting_up_api_key.html)
to know how to apply and use the API key with `rromeo`.

## Usage

`rromeo` contains functions to retrieve data from the SHERPA/RoMEO API
(for a complete overview please refer to the
[vignette](https://docs.ropensci.org/rromeo/articles/rromeo.html)). The
data is released under the [Creative Commons
Attribution-NonCommercial-ShareAlike 2.5 (CC BY-NC-SA 2.5)
license](https://creativecommons.org/licenses/by-nc-sa/2.5/). A
suggestion of citation is included in `rromeo` via `citation("rromeo")`.

`rromeo` functions are prefixed with `rr_` such as `rr_journal_name()`
that lets you retrieve a journal policy information using the title of a
journal:

``` r
rromeo::rr_journal_name("Journal of Biogeography", qtype = "exact")
#>                     title provided_issn      issn romeocolour preprint
#> 1 Journal of Biogeography          <NA> 0305-0270      yellow      can
#>    postprint    pdf pre_embargo post_embargo pdf_embargo
#> 1 restricted cannot        <NA>    12 months        <NA>
```

the `qtype` argument indicates the type of query to make (`exact` for
exact matching of the title, `contains` for partial matching and `starts
with` to match only the beginning of the title).

You can also retrieve a journal information using its ISSN:

``` r
rromeo::rr_journal_issn("0305-0270")
#>                     title provided_issn      issn romeocolour preprint
#> 1 Journal of Biogeography     0305-0270 0305-0270      yellow      can
#>    postprint    pdf pre_embargo post_embargo pdf_embargo
#> 1 restricted cannot        <NA>    12 months        <NA>
```

`rromeo` also provides a function to retrieve information based on
publisher ID `rr_publisher()`.

SHERPA/RoMEO provides a synthetic “colour” for each journal, the colour
summarizes the editorial policy of a journal:

| RoMEO colour | Archiving policy                                        |
| :----------- | :------------------------------------------------------ |
| `green`      | can archive preprint, postprint and publisher’s version |
| `blue`       | can archive postprint **or** publisher’s version        |
| `yellow`     | can archive preprint                                    |
| `white`      | archiving not formally supported                        |

(Table taken from
<http://www.sherpa.ac.uk/romeo/definitions.php#colours>)

`rromeo` lets you retrieve the policies of all journals of a given
colour using the function `rr_romeo_colour()` (**NOTE:** this function
can be slow as there are many journals to retrieve):

``` r
green_journals = rromeo::rr_romeo_colour("green")
green_journals[8:12,]
#>    romeoid                                                   publisher
#> 8     1128 Association for Information Science and Technology (ASIS&T)
#> 9     1937                                       University of Arizona
#> 10    2951                               Geological Society of America
#> 11    2521                              University of California Press
#> 12    2306                                  Optical Society of America
#>                  alias romeocolour preprint postprint        pdf
#> 8              JASIS&T       green      can       can     cannot
#> 9          Radiocarbon       green      can       can restricted
#> 10           GSA Today       green      can       can        can
#> 11            Collabra       green      can       can        can
#> 12 No Paid Open Access       green      can       can     cannot
```

## Dependency network (Imports only)

<img src="man/figures/README-dependency_network_imports-1.svg" width="100%" />

## Dependency network (Imports and Suggests)

<img src="man/figures/README-dependency_network_full-1.svg" width="100%" />

## Contributing to `rromeo`

We welcome contribution to `rromeo`\! Please read the [contribution
guidelines](https://docs.ropensci.org/rromeo/CONTRIBUTING.html) if you
want to contribute, as well as the below-mentioned Code of Conduct.

## Code of Conduct

Please note that the `rromeo` project is released with a [Contributor
Code of Conduct](https://ropensci.org/code-of-conduct/).
By contributing to this project, you agree to abide by its terms.

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# rromeo 0.1.1

## NEW FEATURES

* add column for provided ISSNs in `rr_journal_*()` functions (fix #54)

## BUG FIXES

* Move to static vignettes to avoid issues with CRAN and others
* Update cassette following `vcr` 0.5.0 update

# rromeo 0.1.0

## NEW FEATURES

* released to CRAN
## Test environments
* local Windows 8.1, R-release
* Ubuntu 16.04 (on travis-ci), R-oldrel
* Ubuntu 16.04 (on travis-ci), R-release
* Ubuntu 16.04 (on travis-ci), R-devel
* Windows Server 2012 R2, R-release (AppVeyor CI)
* win-builder R-oldrel
* win-builder R-release
* win-builder R-devel
* Ubuntu Linux 16.04 LTS, R-release, GCC (rhub)
* Fedora Linux, R-devel, clang, gfortran (rhub)
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit (rhub)
* Debian Linux, R-patched, GCC (rhub)

## R CMD check results

0 errors | 0 warnings | 0 note
<!-- Please use a feature branch (i.e., put your work in a new branch that has a name that reflects the feature you are working on; https://docs.gitlab.com/ee/workflow/workflow.html) -->

<!-- If authentication is involved: do not share your username/password, or api keys/tokens in this pull request - most likely the maintainer will have their own equivalent key -->

<!-- If you've updated a file in the man-roxygen directory, make sure to update the man/ files by running devtools::document() or similar as .Rd files should be affected by your change -->

<!--- Provide a general summary of your changes in the Title above -->

## Description
<!--- Describe your changes in detail -->

## Related Issue
<!--- if this closes an issue make sure include e.g., "fix #4"
or similar - or if just relates to an issue make sure to mention
it like "#4" -->

## Example
<!--- if introducing a new feature or changing behavior of existing
methods/functions, include an example if possible to do in brief form -->

<!--- Did you remember to include tests? Unless you're just changing
grammar, please include new tests for your change -->

# Contributing to `rromeo`

We welcome contributions to `rromeo`!
This outlines how to propose a change to rromeo.

### Fixing typos

Small typos or grammatical errors in documentation may be edited directly using
the GitHub web interface, so long as the changes are made in the _source_ file.

*  YES: you edit a roxygen comment in a `.R` file below `R/`.
*  NO: you edit an `.Rd` file below `man/`.

### Prerequisites

Before you make a substantial pull request, you should always file an issue and
make sure someone from the team agrees that it’s a problem. If you’ve found a
bug, create an associated issue and illustrate the bug with a minimal 
[reprex](https://www.tidyverse.org/help/#reprex).

### Pull request process

*  We recommend that you create a Git branch for each pull request (PR).  
*  Look at the Travis and AppVeyor build status before and after making changes.
The `README` should contain badges for any continuous integration services used
by the package.  
*  We use [roxygen2](https://cran.r-project.org/package=roxygen2), with
[Markdown syntax](https://cran.r-project.org/web/packages/roxygen2/vignettes/markdown.html), 
for documentation.  
*  We use [testthat](https://cran.r-project.org/package=testthat). Contributions
with test cases included are easier to accept.  
*  For user-facing changes, add a bullet to the top of `NEWS.md` below the
current development version header describing the changes made followed by your
GitHub username, and links to relevant issue(s)/PR(s).

### Code of Conduct

Please note that the rromeo project is released with a
[Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By 
contributing to this project you agree to abide by its terms.
<!-- If authentication is involved: do not share your username/password, or api keys/tokens in this issue - most likely the maintainer will have their own equivalent key -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
