brranching
==========



[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-check](https://github.com/ropensci/brranching/workflows/R-check/badge.svg)](https://github.com/ropensci/brranching/actions/)
[![codecov.io](https://codecov.io/github/ropensci/brranching/coverage.svg?branch=master)](https://codecov.io/github/ropensci/brranching?branch=master)
[![cran checks](https://cranchecks.info/badges/worst/brranching)](https://cranchecks.info/pkgs/brranching)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/brranching)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/brranching)](https://cran.r-project.org/package=brranching)

## Description

Brranching is an interface to many different sources of phylogenetic data (currently only from Phylomatic (http://phylodiversity.net/phylomatic/), but more sources to come) that allows users to query for phylogenetic data using taxonomic names.  

For `brranching::phylomatic_names()` function you should get an NCBI Entrez API key. NCBI Entrez doesn't require that 
you use an API key, but you get higher rate limit with a key, from 3 to 10 requests per second, so do 
get one. Run `taxize::use_entrez()` or see https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/
for instructions.

## Installation

Stable CRAN version


```r
install.packages("brranching")
```

Or dev version


```r
remotes::install_github("ropensci/brranching")
```


```r
library("brranching")
```

## Phylomatic


```r
taxa <- c("Poa annua", "Phlox diffusa", "Helianthus annuus")
tree <- phylomatic(taxa=taxa, get = 'POST')
plot(tree, no.margin=TRUE)
```

![plot of chunk unnamed-chunk-5](tools/img/unnamed-chunk-5-1.png)

You can pass in up to about 5000 names. We can use `taxize` to get a random set of plant species names.


```r
library("taxize")
spp <- names_list("species", 200)
out <- phylomatic(taxa = spp, get = "POST")
plot(out, show.tip.label = FALSE)
```

![plot of chunk unnamed-chunk-6](tools/img/unnamed-chunk-6-1.png)

## Bladj


```r
library("phylocomr")
ages_df <- data.frame(
  a = c('malpighiales','eudicots','ericales_to_asterales','plantaginaceae',
        'malvids', 'poales'),
  b = c(81, 20, 56, 76, 47, 71)
)
phylo_file <- system.file("examples/phylo_bladj", package = "phylocomr")
phylo_str <- readLines(phylo_file)
x <- rbladj(tree = phylo_str, ages = ages_df)
library(ape)
plot(x)
```

![plot of chunk unnamed-chunk-7](tools/img/unnamed-chunk-7-1.png)

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/brranching/issues).
* License: MIT
* Get citation information for `brranching` in R doing `citation(package = 'brranching')`
* Please note that this project is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By participating in this project you agree to abide by its terms.
brranching 0.7.0
================

### MINOR IMPROVEMENTS

* vignettes fix (#45)


brranching 0.6.0
================

### MINOR IMPROVEMENTS

* the APG dataset in `taxize` package was updated in taxize `v0.9.97` - changes made to comply with the changes in the APG dataset structure  (#41)


brranching 0.5.0
================

### MINOR IMPROVEMENTS

* now using package `conditionz` in the `phylomatic_names()` function for handling messages from the `taxize` package about the user not having an API key set  (#36) (#40)

### BUG FIXES

* require newest `phylocomr` version that has fixes for various mis-behavior  (#38) (#39)


brranching 0.4.0
================

### NEW FEATURES

* in the `phylomatic_local()` function now using `phylocomr::ph_phylomatic`  instead of shelling out to Phylocom via `system`. A number of parameters are gone due to the change internally (#30) (#35)
* in the `rbladj()` function now using `phylocomr::ph_bladj` instead of shelling out to Phylocom via `system` (#30) (#35)
* added a package vignette (#31) (#34) thanks @fozy81
* added new dataset of four phylogenetic trees that can be used in `phylomatic_local()`, see `?phylomatic_trees`

### MINOR IMPROVEMENTS

* added docs to `phylomatic_names()` and the README on using NCBI Entrez API keys


brranching 0.3.0
================

### NEW FEATURES

* gains new function `bladj`  (#18)
* replaced `httr` with `crul` for http requests (#25)

### MINOR IMPROVEMENTS

* fix links to readme images (#29) (#26)
* `verbose` param in `phylomatic()` function changed to `mssgs`


brranching 0.2.0
================

### NEW FEATURES

* Added function `phylomatic_local()` to use Phylomatic locally. 
Phylomatic is a set of Awk scripts, which have to be downloaded 
by the user. After downloading, this function uses the local version 
of Phylomatic (Same as that that runs as a web service). This is 
advantageous especially when dealing with large queries. (#13)

### MINOR IMPROVEMENTS

* Fixed `clean` parameter in `phylomatic()` and `phylomatic_local()`
to expect a logical (`TRUE` or `FALSE`) instead of a "true" or "false". (#15)
* A related change to that above, changed reading newick strings to use 
`phytools::read.newick()` instead of `ape::read.tree()`, which handles
the result of `clean=FALSE` in `phylomatic()` and `phylomatic_local()` (#16)
* Documented that in the `storedtree` parameter of `phylomatic()` and 
`phylomatic_local()` the tree from Zanne et al. is also available by using
`storedtree="zanne2014"` (#19)


brranching 0.1.0
================

### NEW FEATURES

* Released to CRAN.
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
https://contributor-covenant.org/version/1/0/0/
## Test environments

* local macOS install, R 4.0.5 Patched
* ubuntu 16.04 (on GitHub Actions), R 4.0.5
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are no reverse dependencies.

---

This version fixes the rmarkdown/markdown dependency issue for vignettes that Kurt emailed maintainers about.

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

* Submit an issue on the [Issues page](https://github.com/ropensci/brranching/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/brranching.git`
* Make sure to track progress upstream (i.e., on our version of `brranching` at `ropensci/brranching`) by doing `git remote add upstream https://github.com/ropensci/brranching.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Push up to your account
* Submit a pull request to home base at `ropensci/brranching`

### Also, check out our [discussion forum](https://discuss.ropensci.org)

## Prefer to Email? Get in touch: [myrmecocystus@gmail.com](mailto:myrmecocystus@gmail.com)

### Thanks for contributing!
<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below. If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{brranching vignette}
%\VignetteEncoding{UTF-8}
-->

brranching - an interface to phylogenetic data
======



## Description

Brranching is an interface to many different sources of phylogenetic data (currently only from [Phylomatic](http://phylodiversity.net/phylomatic/), but more sources to come). It is used to query for phylogenetic data using taxonomic names and can be used to visualise the evolutionary history and relationships among individuals or groups of organisms. 

## Installation

Stable CRAN version


```r
install.packages("brranching")
```

Or dev version


```r
install.packages("devtools")
devtools::install_github("ropensci/brranching")
```


```r
library("brranching")
```

## Phylomatic

Query Phylomatic for a phylogenetic tree.


```r
taxa <- c("Poa annua", "Phlox diffusa", "Helianthus annuus")
tree <- phylomatic(taxa=taxa, get = 'POST')
plot(tree, no.margin=TRUE)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png)

You can pass in up to about 5000 names. We can use [`taxize`](https://github.com/ropensci/taxize/) to get a random set of plant species names.


```r
library("taxize")
spp <- names_list("species", 200)
out <- phylomatic(taxa = spp, get = "POST")
plot(out, show.tip.label = FALSE)
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png)

## Bladj


```r
library("phylocomr")
ages_df <- data.frame(
  a = c('malpighiales','eudicots','ericales_to_asterales','plantaginaceae',
        'malvids', 'poales'),
  b = c(81, 20, 56, 76, 47, 71)
)
phylo_file <- system.file("examples/phylo_bladj", package = "phylocomr")
phylo_str <- readLines(phylo_file)
x <- rbladj(tree = phylo_str, ages = ages_df)
library(ape)
plot(x)
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png)
