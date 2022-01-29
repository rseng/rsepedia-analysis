rgnparser
=========



[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/rgnparser)](https://cranchecks.info/pkgs/rgnparser)
[![R-check](https://github.com/ropensci/rgnparser/workflows/R-check/badge.svg)](https://github.com/ropensci/rgnparser/actions/)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/rgnparser)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/rgnparser)](https://cran.r-project.org/package=rgnparser)

**rgnparser**: Parse Scientific Names

Docs: https://ropensci.github.io/rgnparser/

## Installation


```r
install.packages("rgnparser")
# OR
remotes::install_github("ropensci/rgnparser")
```


```r
library("rgnparser")
```

## Install gnparser

The command line tool written in Go, gnparser, is required to use this package.

If you want to install gnparser on your own, instructions can be found at the
gnparser repo (https://github.com/gnames/gnparser)

There is a helper function in **rgnparser** for downloading and installing
gnparser on major operating systems (macOS, Windows, Linux):


```r
rgnparser::install_gnparser()
```

It installs the latest gnparser version by default, but you can specify which 
version to install. You can also install gnparser outside of R yourself
(see above).


## Meta

* Please [report any issues or bugs](https://github.com/ropensci/rgnparser/issues).
* License: MIT
* Get citation information for `rgnparser` in R doing `citation(package = 'rgnparser')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
rgnparser 0.2.0
===============

### NEW FEATURES

* A new gnparser version (v1) is out. In addition, gnparser has moved from Gitlab to Github; which also required changes to `install_gnparser()` because we download Go binaries from the gnparser source repository (#7)
* As part of new gnparser version, arguments have changed in `gn_parse()` and `gn_parse_tidy()`: `format` has been removed. `batch_size` and `ignore_tags` were added to both functions, while `details` was added to `gn_parse()` only. See docs for details.  (#11)
* gnparser v1 or greater is now required (#10)

### DEFUNCT

* `gn_debug()` is now defunct. the gnparser command for this function was removed in gnparser v1 (#9)

### BUG FIXES

* `gn_version()` was broken with the new gnparser version, fixed now (#8)
* xxx (#xx)


rgnparser 0.1.0
===============

### NEW FEATURES

* First submission to CRAN.
## Test environments
* local R installation, R 4.0.3
* ubuntu 16.04 (on GitHub Actions), R 4.0.3
* win-builder (release and devel)

## R CMD check results

0 errors | 0 warnings | 0 notes

-----

This submission makes appropriate changes for a new major version of the gnparser Go command line tool.

Thanks!
Scott Chamberlain
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

# CONTRIBUTING #

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
*  We recommend the tidyverse [style guide](http://style.tidyverse.org).
You can use the [styler](https://CRAN.R-project.org/package=styler) package to
apply these styles, but please don't restyle code that has nothing to do with 
your PR.  
*  We use [roxygen2](https://cran.r-project.org/package=roxygen2).  
*  We use [testthat](https://cran.r-project.org/package=testthat). Contributions
with test cases included are easier to accept.  
*  For user-facing changes, add a bullet to the top of `NEWS.md` below the
current development version header describing the changes made followed by your
GitHub username, and links to relevant issue(s)/PR(s).

### Code of Conduct

Please note that the `pegax` project is released with a
[Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project you agree to abide by its terms.

### See rOpenSci [contributing guide](https://devguide.ropensci.org/contributingguide.html)
for further details.

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.
<!-- If authentication is involved: do not share your username/password, or api keys/tokens in this issue - most likely the maintainer will have their own equivalent key -->

<!-- Do not share screen shots of code. Share actual code in text format. -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/). If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
