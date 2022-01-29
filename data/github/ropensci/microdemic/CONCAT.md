microdemic
==========



[![cran checks](https://cranchecks.info/badges/worst/microdemic)](https://cranchecks.info/pkgs/microdemic)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-check](https://github.com/ropensci/microdemic/workflows/R-check/badge.svg)](https://github.com/ropensci/microdemic/actions?query=workflow%3AR-check)
[![codecov](https://codecov.io/gh/ropensci/microdemic/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/microdemic)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/microdemic)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/microdemic)](https://cran.r-project.org/package=microdemic)

`microdemic` - Microsoft Academic Client

Web interface: https://academic.microsoft.com/

API docs:
- https://docs.microsoft.com/en-us/azure/cognitive-services/academic-knowledge/
- https://msr-apis.portal.azure-api.net/docs/services/academic-search-api/operations/565d9001ca73072048922d97

Get a API key at https://msr-apis.portal.azure-api.net/signin

## install

cran version


```r
install.packages("microdemic")
```

dev version


```r
remotes::install_github("ropensci/microdemic")
```


```r
library("microdemic")
```

## Evaluate

See the [query expression syntax](https://docs.microsoft.com/en-us/azure/cognitive-services/academic-knowledge/queryexpressionsyntax)
for help on how to construct queries - for this and other functions


```r
ma_evaluate(query = "Y='19'...")
#> # A tibble: 10 x 8
#>    logprob     prob      Id Ti                         Y    CC AA      J.JN     
#>      <dbl>    <dbl>   <dbl> <chr>                  <int> <int> <list>  <chr>    
#>  1   -13.7  1.10e-6  2.12e9 image forming device    1992 32430 <df[,1… <NA>     
#>  2   -13.8  1.04e-6  1.86e9 standard methods for …  1992 81915 <df[,1… <NA>     
#>  3   -13.8  1.03e-6  2.16e9 the nature of statist…  1995 30098 <df[,1… <NA>     
#>  4   -13.9  9.18e-7  2.91e9 fuzzy sets              1996 44493 <df[,1… <NA>     
#>  5   -13.9  9.12e-7  2.16e9 gapped blast and psi …  1997 61351 <df[,1… nucleic …
#>  6   -13.9  8.77e-7  2.23e9 manufacture of semico…  1992 29044 <df[,1… <NA>     
#>  7   -14.1  7.44e-7  2.15e9 statistical learning …  1998 21495 <df[,1… <NA>     
#>  8   -14.1  7.22e-7  2.12e9 neural networks a com…  1998 24498 <df[,1… <NA>     
#>  9   -14.2  6.96e-7  1.98e9 generalized gradient …  1996 84892 <df[,1… physical…
#> 10   -14.2  6.86e-7  2.99e9 particle swarm optimi…  1995 12985 <df[,1… <NA>
```

## Calchistogram


```r
res <- ma_calchist(query = "And(Composite(AA.AuN=='jaime teevan'),Y>2012)",
   atts = c('Y', 'F.FN'))
res$histograms$histogram
#> [[1]]
#>   value   logprob count
#> 1  2013 -17.01346    19
#> 2  2014 -17.07550    14
#> 3  2015 -17.42947    15
#> 4  2016 -17.50792    17
#> 5  2019 -17.65841     6
#> 6  2017 -18.09004    11
#> 7  2018 -18.47003     7
#> 8  2020 -18.95307     4
#> 
#> [[2]]
#>                         value   logprob count
#> 1            computer science -15.66718    75
#> 2              world wide web -16.54773    28
#> 3               crowdsourcing -16.63486    24
#> 4  human computer interaction -16.85747    19
#> 5               search engine -17.00135    14
#> 6       information retrieval -17.29366    11
#> 7                  multimedia -17.79757    10
#> 8     artificial intelligence -17.81949     6
#> 9            search analytics -17.83150     5
#> 10               data science -17.90774    11
```

## Abstract


```r
ma_abstract(query = "Y='19'...", count = 5)
#> # A tibble: 4 x 2
#>           Id abstract                                                           
#>        <dbl> <chr>                                                              
#> 1 2119113870 An image forming device has: an image forming body on which an ima…
#> 2 1856219842 Set your standards with these standard methods. This is it: the mo…
#> 3 2156909104 Setting of the learning problem consistency of learning processes …
#> 4 2158714788 The BLAST programs are widely used tools for searching protein and…
```


## Meta

* Please [report any issues or bugs](https://github.com/ropensci/microdemic/issues).
* License: MIT
* Get citation information for `microdemic` in R doing `citation(package = 'microdemic')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
microdemic 0.6.0
================

### NEW FEATURES

* The `E` attribute we used heavily for data for various functions here has been dropped and its sub fields brought to the top level. Just internal changes - except for that the default fields in `ma_calchist`, `ma_evaluate`, and `ma_search` now do not include `E`. See the github issue for link for details  (#17)

### DEFUNCT

* `ma_similarity()` is now defunct as the service at least appears to be gone or down often enough its not worth supporting (#16)
* `ma_graph_search()` is now defunct as the service at least appears to be gone or down often enough its not worth supporting (#18)

microdemic 0.5.0
================

### NEW FEATURES

* new author Christopher Baker (#13)
* docs website added

### MINOR IMPROVEMENTS

* all functions now pass along detailed error messages; before we were erroring but just giving the http status code and generic message; as part of this fix now importing httpcode pkg (#7) (#12) (#14)
* use fake API key on travis to avoid test suite error (#11)
* internals of `ma_search()` have changed a bit; get in touch if you have any questions about this function

### BUG FIXES

* fix to internal function `invabs2abs()` within `ma_abstract()`: the inverted index returned from Microsoft can have missing values and we were incorrectly inserting NA's into those spots resulting in NA's in abstracts (#8) (#9)

microdemic 0.4.0
================

### MINOR IMPROVEMENTS

* do testing with `vcr` (#5) (#6)
* make the output of `ma_abstract()` a data.frame and add `Id` column to 
facilitate matching to other data (#4)

microdemic 0.3.0
================

### MINOR IMPROVEMENTS

* Improve docs on how to use and get API keys (#3)
* Change base URL from `westus.api.cognitive.microsoft.com` to `api.labs.cognitive.microsoft.com`. Because of this change, you need to get use an API key from a the microsoft labs website. Get a key at some url and see `?microdemic-package` for details on how to use it (#2)


microdemic 0.2.0
================

### NEW FEATURES

* Many of the functions gain a new parameter `model` with value of 
'latest' or 'beta-2015'. 


microdemic 0.1.0
================

### NEW FEATURES

* released to CRAN
## Test environments

* local OS X install, R 4.0.3 Patched
* ubuntu 16.04 (on travis-ci), R 4.0.3
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

* I have run R CMD check on the 1 reverse dependency. No problems were
found. Summary at <https://github.com/ropensci/microdemic/tree/master/revdep>

------

This version makes two functions defunct and makes changes aligning with remote API changes.

Thanks!
Scott Chamberlain
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

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

### Bugs?

* Submit an issue on the [Issues page](https://github.com/ropensci/microdemic/issues)

### Issues and Pull Requests

If you are considering a pull request, you may want to open an issue first to discuss with the maintainer(s).

### Code contributions

* Fork this repo to your GitHub account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/microdemic.git`
* Make sure to track progress upstream (i.e., on our version of `microdemic` at `ropensci/microdemic`) by doing `git remote add upstream https://github.com/ropensci/microdemic.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch - see <https://guides.github.com/introduction/flow/> for how to contribute by branching, making changes, then submitting a pull request)
* Push up to your account
* Submit a pull request to home base (likely master branch, but check to make sure) at `ropensci/microdemic`

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.

### Prefer to Email? 

Do not email, open an issue.
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```

</details>
# Platform

|field    |value                                       |
|:--------|:-------------------------------------------|
|version  |R version 4.0.3 Patched (2020-11-02 r79396) |
|os       |macOS Catalina 10.15.7                      |
|system   |x86_64, darwin17.0                          |
|ui       |X11                                         |
|language |(EN)                                        |
|collate  |en_US.UTF-8                                 |
|ctype    |en_US.UTF-8                                 |
|tz       |US/Pacific                                  |
|date     |2020-11-19                                  |

# Dependencies

|package    |old   |new        |Δ  |
|:----------|:-----|:----------|:--|
|microdemic |0.5.0 |0.6.0      |*  |
|crayon     |NA    |1.3.4.9000 |*  |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*