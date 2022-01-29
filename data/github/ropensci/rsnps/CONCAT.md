
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rsnps

[![cran
checks](https://cranchecks.info/badges/worst/rsnps)](https://cranchecks.info/pkgs/rsnps/)
[![R build
status](https://github.com/ropensci/rsnps/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/rsnps/actions)
[![Build
status](https://ci.appveyor.com/api/projects/status/d2lv98726u6t9ut5/branch/master)](https://ci.appveyor.com/project/sckott/rsnps/branch/master/)
[![codecov.io](https://codecov.io/github/ropensci/rsnps/coverage.svg?branch=master)](https://codecov.io/github/ropensci/rsnps?branch=master)
[![cran
version](https://www.r-pkg.org/badges/version/rsnps)](https://cran.r-project.org/package=rsnps)
[![rstudio mirror
downloads](https://cranlogs.r-pkg.org/badges/rsnps?color=E664A4)](https://github.com/r-hub/cranlogs.app)

This package gives you access to data from OpenSNP and NCBI’s dbSNP SNP
database.

## NOTE

`rsnps` used to be `ropensnp`

## Data sources

This set of functions/package accesses data from:

-   openSNP.org
    -   <https://opensnp.org>
    -   See documentation on the openSNP API
        <https://opensnp.org/faq#api>
    -   See blog post about their API
        <https://opensnp.wordpress.com/2012/01/18/some-progress-on-the-api-json-endpoints/>
    -   Relevant functions:
        -   `allgensnp()`, `allphenotypes()`, `annotations()`,
            `download_users()`, `fetch_genotypes()`, `genotypes()`,
            `phenotypes()`, `phenotypes_byid()`, `users()`
-   NCBI’s dbSNP SNP database
    -   See <https://www.ncbi.nlm.nih.gov/snp/> for more details
    -   Relevant function:
        -   `ncbi_snp_query()`

## Install

Install from CRAN

``` r
install.packages("rsnps")
```

Or dev version

``` r
install.packages("remotes")
remotes::install_github("ropensci/rsnps")
```

``` r
library("rsnps")
```

## Usage

### NCBI dbSNP data

``` r
snps <- c("rs332", "rs420358", "rs1837253", "rs1209415715", "rs111068718")
ncbi_snp_query(snps)
```

    #> # A tibble: 4 × 16
    #>   query        chromosome        bp class rsid   gene   alleles ancestral_allele
    #>   <chr>        <chr>          <dbl> <chr> <chr>  <chr>  <chr>   <chr>           
    #> 1 rs332        7          117559593 del   rs121… "CFTR… TTT, d… TTT             
    #> 2 rs420358     1           40341239 snv   rs420… ""     A,C,G,T A               
    #> 3 rs1837253    5          111066174 snv   rs183… ""     T,C     T               
    #> 4 rs1209415715 9           41782316 snv   rs120… ""     T,A,C   T               
    #> # … with 8 more variables: variation_allele <chr>, seqname <chr>, hgvs <chr>,
    #> #   assembly <chr>, ref_seq <chr>, minor <chr>, maf <dbl>,
    #> #   maf_population <list>

### openSNP data

`genotypes()` function

``` r
genotypes('rs9939609', userid='1,6,8', df=TRUE)
```

    #>    snp_name snp_chromosome snp_position                 user_name user_id
    #> 1 rs9939609             16     53786615 Bastian Greshake Tzovaras       1
    #> 2 rs9939609             16     53786615              Nash Parovoz       6
    #> 3 rs9939609             16     53786615         Samantha B. Clark       8
    #>   genotype_id genotype
    #> 1           9       AT
    #> 2           5       AT
    #> 3           2       TT

`phenotypes()` function

``` r
out <- phenotypes(userid=1)
out$phenotypes$`Hair Type`
```

    #> $phenotype_id
    #> [1] 16
    #> 
    #> $variation
    #> [1] "straight"

For more detail, see the [vignette: rsnps
tutorial](https://github.com/ropensci/rsnps/tree/master/vignettes).

## Meta

-   Please [report any issues or
    bugs](https://github.com/ropensci/rsnps/issues/).
-   License: MIT
-   Get citation information for `rsnsps` in R doing
    `citation(package = 'rsnps')`
-   Please note that this package is released with a [Contributor Code
    of Conduct](https://ropensci.org/code-of-conduct/). By contributing
    to this project, you agree to abide by its terms.

[![ropensci\_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
rsnps 0.5.0
===========

### NEW FEATURES

* ncbi_snp_query(): enable allele frequency for different reference populations, ncbi_snp_query() outputs now a tibble (#97).

### MINOR IMPROVEMENTS

* ncbi_snp_query(): replace calls to RJSONIO with equivalent in jsonlite (#98).
* unit tests for ncbi_snp_query(): added a tolerance to any allele frequency checks.
* move vignette source to /vignettes and precompute using an R script.

### DOCUMENTATION FIXES

* Updated vignette.

rsnps 0.4.0
===========

### MAJOR IMPROVEMENTS

NCBI / dbSNP changed their API:

* Rewrote `ncbi_snp_query` to accommodate the new API (#86, #88). 
* Removed the functions `ncbi_snp_query2` an `ncbi_snp_summary`. 

### MINOR IMPROVEMENTS

* Reordered `ncbi_snp_query` dataframe output to have chromosome and bp beside each other (#70).
* Changed `ncbi_snp_query` parameter (`SNPs`) to lower case (`snps`). 

### DOCUMENTATION FIXES

* Restructured and fixed a typo in `README.Rmd` and added link to vignette (#63).
* Added info of two new maintainers to `DESCRIPTION`. 
* Added relevant API links to vignette. 

### BUG FIXES

* Fixed the test for `allphenotypes` function by making it less specific (#72). 


rsnps 0.3.0
===========

### DEPRECATED AND DEFUNCT

* `ld_search()` is now defunct. The Broad Institute has taken down the SNAP service behind the function. (#46) (#53) (#60)

### NEW FEATURES

* the three NCBI functions gain a new parameter `key` for passing in an NCBI Entrez API key. You can alternatively (and we encourage this) store your key as an environment variable and we'll use that instead. The key allows you to have higher rate limits than without a key (#58)
* gains new function `ncbi_snp_summary()` for summary data on a SNP (#31)

### MINOR IMPROVEMENTS

* http requests are now done using `crul` instead of `httr` (#44)
* now using markdown formatted documentation (#56)
* documented in `ncbi_snp_query()` that we can not change the assembly (#49)

### BUG FIXES

* fix to `ncbi_snp_query2()`: when many IDs passed in, we were failing with a "URI too long" message. We now check how many Ids are passed in and do a POST request as needed  (#39)
* fixed problem in `ncbi_snp_query()` where it wasn't pulling out correctly the gene name and BP position (#25)



rsnps 0.2.0
===========

### NEW FEATURES

* `LDSearch()` is now `ld_search()`, but `LDSearch()` still works until 
the next CRAN release when it will be defunct (#33)
* `NCBI_snp_query()` is now `ncbi_snp_query()`, but `NCBI_snp_query()` still 
works until the next CRAN release when it will be defunct (#33)
* `NCBI_snp_query2()` is now `ncbi_snp_query2()`, but `NCBI_snp_query2()` still 
works until the next CRAN release when it will be defunct (#33)

### MINOR IMPROVEMENTS

* Namespace all base R package function calls (#21)
* Improve `httr::content` call to parse to text, and `encoding = "UTF-8"` 
(#24)
* Added tests for `ld_search()` (#12)
* Added tests for `ncbi_snp_query()` and `ncbi_snp_query2()` (#13)
* Added ancestral allele output to `ncbi_snp_query()` (#23)

### BUG FIXES

* Fix to `fetch_genotypes()`, was failing sometimes when the commented
metadata lines at top varied in length (#22)
* Fix to `ld_search()` (#32)


rsnps 0.1.6
===========

### MINOR IMPROVEMENTS

* All examples now in `\dontrun`. (#11)
* Added additional tests for `LDSearch()` and `NCBI_snp_query()`.
* Added a vignette.

### BUG FIXES

* Bugs fixed in `LDSearch()`, which were actually bugs in `NCBI_snp_query()`. (#9)
* Bug fixed in `NCBI_snp_query()` as chromosome might also be "X". 

rsnps 0.1.0
===========

### NEW FEATURES 

* Bug fixes to all openSNP functions.

rsnps 0.0.5
===========

### NEW FEATURES 

* released to CRAN
## Test environments

* local OS X install, R 4.1.1
* Ubuntu Linux 20.04.1 LTS (on R-hub), R 4.1.2
* Fedora Linux (on R-hub) R-devel
* Windows (devel and release)

## R CMD check results

There were no ERRORs or WARNINGs. 

There is one NOTE that is only found on Windows (Server 2022, R-devel 64-bit): 

```
* checking for detritus in the temp directory ... NOTE
Found the following files/directories:
  'lastMiKTeXException'
```
As noted in [R-hub issue #503](https://github.com/r-hub/rhub/issues/503), this could be due to a bug/crash in MiKTeX and can likely be ignored.

## Reverse dependencies

* We have run R CMD check on the 1 downstream dependency
(<https://github.com/ropensci/rsnps/blob/master/revdep/README.md>).
No problems were found related to this package.

---

This version includes a new feature and two minor improvement of the function `ncbi_snp_query`,
and a change in the vignette (we now pre-compile the vignette to avoid long runtimes). 


Thanks!
Julia Gustavsen and Sina Rüeger
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

### Please contribute!

We love collaboration.

### Bugs?

* Submit an issue on the [Issues page](https://github.com/{owner}/{repo}/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/{repo}.git`
* Make sure to track progress upstream (i.e., on our version of `{repo}` at `{owner}/{repo}`) by doing `git remote add upstream https://github.com/{owner}/{repo}.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Push up to your account
* Submit a pull request to home base (likely master branch, but check to make sure) at `{owner}/{repo}`

### Also, check out our [discussion forum](https://discuss.ropensci.org)

### Prefer to Email? Email the person listed as maintainer in the `DESCRIPTION` file of this repo.

### Thanks for contributing!
---
name: Bug report
about: Create a report to help us improve rsnps
title: ''
labels: 'Bug'
assignees: ''

---

**Please do not submit bug reports for issues with the OpenSNP or NCBI’s dbSNP server itself. rsnps is only an API client, it does not have control over the server hosting the data.**

**Describe the bug**
A clear and concise description of what the bug is.

**Code to reproduce**
A clear set of code illustrating the steps to reproduce the behaviour.

**Expected behaviour**
A clear and concise description of what you expected to happen.

**OS and R versions (please complete the following information):**
 - OS: [e.g., Windows10, macOS, Linux]
 - R Version [e.g., 3.6.2]

**Session Info: `devtools::session_info()` or `sessionInfo()`**

<details> <summary><strong>Session Info</strong></summary>
  
```r

```

**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for rsnps
title: ''
labels: 'Feature'
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. 

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
# Platform

|field    |value                        |
|:--------|:----------------------------|
|version  |R version 3.6.2 (2019-12-12) |
|os       |macOS Mojave 10.14.6         |
|system   |x86_64, darwin15.6.0         |
|ui       |RStudio                      |
|language |(EN)                         |
|collate  |en_CA.UTF-8                  |
|ctype    |en_CA.UTF-8                  |
|tz       |Europe/Zurich                |
|date     |2020-06-21                   |

# Dependencies

|package |old   |new        |Δ  |
|:-------|:-----|:----------|:--|
|rsnps   |0.3.0 |0.3.2.9121 |*  |
|plyr    |NA    |1.8.6      |*  |
|xml2    |NA    |1.3.2      |*  |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*