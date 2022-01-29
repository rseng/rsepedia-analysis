rbace
=====



[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/rbace)](https://cranchecks.info/pkgs/rbace)
[![R-check](https://github.com/ropensci/rbace/workflows/R-check/badge.svg)](https://github.com/ropensci/rbace/actions?query=workflow%3AR-check)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/rbace?color=C9A115)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/rbace)](https://cran.r-project.org/package=rbace)


Client for interacting with the Bielefeld Academic Search Engine API.

Docs: https://docs.ropensci.org/rbace/

BASE API docs: https://www.base-search.net/about/download/base_interface.pdf

Access: The BASE API is IP address AND user-agent (see note below) restricted. The user agent is set correctly if you use this package, but you still need to get your IP address(es) white-listed by BASE. Request access at: https://www.base-search.net/about/en/contact.php - Note: the BASE website has a search portal you can use from anywhere; it's just the API that is IP and user-agent restricted.

Terminology:

- an IP address is the numeric label identifying a computer or server. the IP address for a computer can change, e.g., if you connect to a VPN
- a user-agent is a string of text that identifies the software requesting data from a server (in this case BASE's API).

Data from BASE (Bielefeld Academic Search Engine) https://www.base-search.net

[<img src="man/figures/BASE_search_engine_logo.svg.png" width="300">](https://www.base-search.net)

## Install


```r
install.packages("rbace")
```

or the dev version


```r
remotes::install_github("ropensci/rbace")
# OR the below should install the same thing
install.packages("rbace", repos = "https://dev.ropensci.org")
```


```r
library("rbace")
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/rbace/issues).
* License: MIT
* Get citation information for `rbace` in R doing `citation(package = 'rbace')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
rbace 0.2.2
===========

### MINOR IMPROVEMENTS

* improved some tests


rbace 0.2.0
===========

### MINOR IMPROVEMENTS

First version to CRAN.

* replace `dplyr` with `data.table::rbindlist` (#6)
* two new functions: `bs_profile()`, `bs_repositories()` (#9)
* add faceting to `bs_search()` (#10)
* `bs_search()` now uses RETRY/GET http requests (#14)
* `bs_search()` gains new parameter `filter` that's passed to `fq` param for the Solr server (#15)
* `filter` parameter in `bs_search()` can be `character` and `AsIs` in case you need to prevent html escaping (#16)
## Test environments
* local R installation, R 4.0.3
* ubuntu 16.04 (on travis-ci), R 4.0.3
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

---

This version improves some tests.

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

* Submit an issue on the [Issues page](https://github.com/ropensci/rbace/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/rbace.git`
* Make sure to track progress upstream (i.e., on our version of `rbace` at `ropensci/rbace`) by doing `git remote add upstream https://github.com/ropensci/rbace.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Push up to your account
* Submit a pull request to home base at `ropensci/rbace`

### Also, check out our [discussion forum](https://discuss.ropensci.org)
<!-- If authentication is involved: do not share your username/password, or api keys/tokens in this issue - most likely the maintainer will have their own equivalent key -->

<!-- Do not share screen shots of code. Share actual code in text format. -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/). If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
