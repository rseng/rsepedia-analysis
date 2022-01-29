rdryad
======



[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-check](https://github.com/ropensci/rdryad/workflows/R-check/badge.svg)](https://github.com/ropensci/rdryad/actions?query=workflow%3AR-check)
[![codecov](https://codecov.io/gh/ropensci/rdryad/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/rdryad)
[![cran checks](https://cranchecks.info/badges/worst/rdryad)](https://cranchecks.info/pkgs/rdryad)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/rdryad)](https://github.com/metacran/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/rdryad)](https://cran.r-project.org/package=rdryad)

`rdryad` is a package to interface with the Dryad data repository.

General Dryad API documentation: https://datadryad.org/api/v2/docs/

rdryad docs: https://docs.ropensci.org/rdryad/

## Installation

Install Dryad from CRAN


```r
install.packages("rdryad")
```

development version:


```r
remotes::install_github("ropensci/rdryad")
```


```r
library('rdryad')
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/rdryad/issues).
* License: MIT
* Get citation information for `rdryad` in R doing `citation(package = 'rdryad')`
* Please note that this project is released with a [Contributor Code of Conduct][coc]. By participating in this project you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)

### Data provided by...

Data is provided from the Dryad API.

[coc]: https://github.com/ropensci/rdryad/blob/master/CODE_OF_CONDUCT.md
rdryad 1.0.0
============

### BREAKING CHANGES

* Package redone to work with new Dryad API v2. Most functions are defunct, and there's three sets of new functions following the three major sets of API routes for datasets, versions, and files. See `?rdryad` for more (#28) (#29)


rdryad 0.4.0
============

### NEW FEATURES

* gains new function `dryad_metadata()` to download Dryad file metadata 
* gains new function `dryad_package_dois()` to get file DOIs for a Dryad package DOI (a package can have many files) (#22)

### MINOR IMPROVEMENTS

* `dryad_files` (formerly `download_url()`) now scrapes Dryad page to get URLs to Dryad files instead of using their API, which was not dependable (#26)
* `dryad_fetch` gains a parameter `try_file_names` (a boolean) which if `TRUE` we try to extract file names out of URLs (#26)

### BUG FIXES

* fix to solr `rdryad` functions to hard code use of `xml` return format, and followlocation to follow any redirects (#27)

### DEFUNCT

* `download_url()` is now defunct, see `dryad_files()`

### NOTE

* two new pacakage dependencies: `tibble` and `data.table`

rdryad 0.3.0
============

### NEW FEATURES

* Move to using `solrium` package instead of `solr` package
for interaction with Dryad's Solr backend (#21) (#24)
* Now using `crul` instead of `httr` for HTTP requests (#23)
* gains two new functions `handle2doi` and `doi2handle` to
convert between handles and DOIs, and DOIs and handles,
respectively (#25)
* `download_url` function name has been changed to `dryad_files`, but
you can still use `download_url` until the next version. In addition,
`download_url`/`dryad_files` parameters `id` is changed to `doi`.

### MINOR IMPROVEMENTS

* `dryad_fetch` is improved, and uses `curl::curl_download` instead of
`download.file`. It now accepts >1 input URL, but `destile` length must
equal number of urls.


rdryad 0.2.0
============

### NEW FEATURES

* Re-worked most of the package.
* New package API, some methods are the same, but many are different. (#16)
* New functions (see functions starting with `d_*()`) to interact
with Dryad Solr search engine (#10)
* OAI-PMH functions now using internally the `oai` package. (#14)

### MINOR IMPROVEMENTS

* Slimmed down dependencies to a smaller set.
* Changed license from CC0 to MIT (#17)
* Added more tests (#18)
* Changed function to get files to only download them, and not attempt to
read them into R, which introduces a very long dependency chain (#15)


rdryad 0.1.1
============

### BUG FIXES

* removed read.jpeg as a dependency


rdryad 0.1
==========

### NEW FEATURES

* released to CRAN
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
(https://contributor-covenant.org), version 1.0.0, available at 
http://contributor-covenant.org/version/1/0/0/
## Test environments

* local OS X install, R 4.0.2
* ubuntu 16.04 (on travis-ci), R 4.0.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are no reverse dependencies.

---

This version completely overhauls the package, most old functions defunct, bump major version for the breaking changes.

Thanks!
Scott Chamberlain
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

### Bugs?

* Submit an issue on the [Issues page](https://github.com/ropensci/rdryad/issues)

### Issues and Pull Requests

If you are considering a pull request, you may want to open an issue first to discuss with the maintainer(s).

### Code contributions

* Fork this repo to your GitHub account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/rdryad.git`
* Make sure to track progress upstream (i.e., on our version of `rdryad` at `ropensci/rdryad`) by doing `git remote add upstream https://github.com/ropensci/rdryad.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch - see <https://guides.github.com/introduction/flow/> for how to contribute by branching, making changes, then submitting a pull request)
* Push up to your account
* Submit a pull request to home base (likely master branch, but check to make sure) at `ropensci/rdryad`

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.
<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
