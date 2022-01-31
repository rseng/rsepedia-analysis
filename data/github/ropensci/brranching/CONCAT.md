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
brranching
==========

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE,
  fig.path="tools/img/"
)
```

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

```{r eval=FALSE}
install.packages("brranching")
```

Or dev version

```{r eval=FALSE}
remotes::install_github("ropensci/brranching")
```

```{r}
library("brranching")
```

## Phylomatic

```{r}
taxa <- c("Poa annua", "Phlox diffusa", "Helianthus annuus")
tree <- phylomatic(taxa=taxa, get = 'POST')
plot(tree, no.margin=TRUE)
```

You can pass in up to about 5000 names. We can use `taxize` to get a random set of plant species names.

```{r}
library("taxize")
spp <- names_list("species", 200)
out <- phylomatic(taxa = spp, get = "POST")
plot(out, show.tip.label = FALSE)
```

## Bladj

```{r}
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

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/brranching/issues).
* License: MIT
* Get citation information for `brranching` in R doing `citation(package = 'brranching')`
* Please note that this project is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By participating in this project you agree to abide by its terms.
<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{brranching vignette}
%\VignetteEncoding{UTF-8}
-->

brranching - an interface to phylogenetic data
======

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE
)
```

## Description

Brranching is an interface to many different sources of phylogenetic data (currently only from [Phylomatic](http://phylodiversity.net/phylomatic/), but more sources to come). It is used to query for phylogenetic data using taxonomic names and can be used to visualise the evolutionary history and relationships among individuals or groups of organisms. 

## Installation

Stable CRAN version

```{r eval=FALSE}
install.packages("brranching")
```

Or dev version

```{r eval=FALSE}
install.packages("devtools")
devtools::install_github("ropensci/brranching")
```

```{r}
library("brranching")
```

## Phylomatic

Query Phylomatic for a phylogenetic tree.

```{r}
taxa <- c("Poa annua", "Phlox diffusa", "Helianthus annuus")
tree <- phylomatic(taxa=taxa, get = 'POST')
plot(tree, no.margin=TRUE)
```

You can pass in up to about 5000 names. We can use [`taxize`](https://github.com/ropensci/taxize/) to get a random set of plant species names.

```{r}
library("taxize")
spp <- names_list("species", 200)
out <- phylomatic(taxa = spp, get = "POST")
plot(out, show.tip.label = FALSE)
```

## Bladj

```{r}
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
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/brranching-package.R
\docType{package}
\name{brranching-package}
\alias{brranching-package}
\alias{brranching}
\title{brranching}
\description{
Phylogenies from many sources
}
\author{
Scott Chamberlain \href{mailto:myrmecocystus@gmail.com}{myrmecocystus@gmail.com}
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phylomatic.R
\name{phylomatic}
\alias{phylomatic}
\title{Query Phylomatic for a phylogenetic tree.}
\usage{
phylomatic(
  taxa,
  taxnames = TRUE,
  get = "GET",
  informat = "newick",
  method = "phylomatic",
  storedtree = "R20120829",
  treeuri = NULL,
  taxaformat = "slashpath",
  outformat = "newick",
  clean = TRUE,
  db = "apg",
  mssgs = TRUE,
  ...
)
}
\arguments{
\item{taxa}{Phylomatic format input of taxa names.}

\item{taxnames}{If \code{TRUE} (default), we get the family names for you to attach
to your species names to send to Phylomatic API. If \code{FALSE}, you have to
provide the strings in the right format.}

\item{get}{'GET' (default) or 'POST' format for submission to the website.}

\item{informat}{One of newick (default), nexml, or cdaordf. If using a stored tree,
informat should always be newick.}

\item{method}{One of phylomatic (default) or convert}

\item{storedtree}{One of R20120829 (Phylomatic tree R20120829 for plants),
smith2011 (Smith 2011, plants), binindaemonds2007 (Bininda-Emonds 2007,
mammals), or zanne2014 (Zanne et al. 2014, plants). Default: R20120829}

\item{treeuri}{URL for a phylogenetic tree in newick format.}

\item{taxaformat}{Only option is slashpath for now. Leave as is.}

\item{outformat}{One of newick, nexml, or fyt.}

\item{clean}{Return a clean tree or not. Default: \code{TRUE}}

\item{db}{One of "ncbi", "itis", or "apg". Default: apg}

\item{mssgs}{Print messages. Default: \code{TRUE}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
Newick formatted tree as \code{phylo} object or
nexml character string
}
\description{
Query Phylomatic for a phylogenetic tree.
}
\details{
Use the web interface at \url{http://phylodiversity.net/phylomatic/}

If you set \code{taxnames = FALSE}, you need to pass in a character
vector, with each element like this example:
\code{"asteraceae/taraxacum/taraxacum_officinale"}, of the form
\code{"family/genus/genus_specfic epithet"}
}
\examples{
\dontrun{
# Input taxonomic names
taxa <- c("Poa annua", "Phlox diffusa", "Helianthus annuus")
tree <- phylomatic(taxa=taxa, get = 'POST')
plot(tree, no.margin=TRUE)

# Genus names
taxa <- c("Poa", "Phlox", "Helianthus")
tree <- phylomatic(taxa=taxa, storedtree='R20120829', get='POST')
plot(tree, no.margin=TRUE)

# Lots of names
taxa <- c("Poa annua", "Collomia grandiflora", "Lilium lankongense", "Phlox diffusa",
"Iteadaphne caudata", "Gagea sarmentosa", "Helianthus annuus")
tree <- phylomatic(taxa=taxa, get = 'POST')
plot(tree, no.margin=TRUE)

# Don't clean - clean=TRUE is default
(tree <- phylomatic(taxa=taxa, clean = FALSE))
## with clean=FALSE, you can get non-splitting nodes, which you
## need to collpase before plotting
library('ape')
plot(collapse.singles(tree), no.margin=TRUE)

# Output NeXML format
taxa <- c("Gonocarpus leptothecus", "Gonocarpus leptothecus", "Lilium lankongense")
out <- phylomatic(taxa=taxa, get = 'POST', outformat = "nexml")
cat(out)

# Lots of names, note that when you have enough names (number depends on length of individual
# names, so there's no per se rule), you will get an error when using `get='GET'`,
# when that happens use `get='POST'`
library("taxize")
spp <- names_list("species", 500)
# phylomatic(taxa = spp, get = "GET")
(out <- phylomatic(taxa = spp, get = "POST", db = "itis"))
plot(out)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rbladj.R
\name{rbladj}
\alias{rbladj}
\title{Run Phylocom's bladj from R}
\usage{
rbladj(tree, ages)
}
\arguments{
\item{tree}{(character/phylo) One of: phylogeny as a newick string (will be
written to a temp file) - OR path to file with a newick
string - OR a an \pkg{ape} \code{phylo} object. required.}

\item{ages}{(data.frame/character) ages data.frame, or path to an ages
file. required.}
}
\value{
Newick formatted tree as \code{phylo} object
}
\description{
Run Phylocom's bladj from R
}
\details{
uses \code{\link[phylocomr:ph_bladj]{phylocomr::ph_bladj()}} under the hood
}
\examples{
\dontrun{
library("phylocomr")

# make an ages data.frame
ages_df <- data.frame(
  a = c('malpighiales','eudicots','ericales_to_asterales','plantaginaceae',
        'malvids', 'poales'),
  b = c(81, 20, 56, 76, 47, 71)
)

# read phylogeny file as a string
phylo_file <- system.file("examples/phylo_bladj", package = "phylocomr")
phylo_str <- readLines(phylo_file)

# Run Bladj, returns phylo object
(x <- rbladj(tree = phylo_str, ages = ages_df))

# load ape and plot tree
library(ape)
plot(x)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/brranching-package.R
\docType{data}
\name{phylomatic_trees}
\alias{phylomatic_trees}
\title{Phylogenies to use with phylomatic}
\format{
A list with 4 character strings:
\itemize{
\item R20120829 - 2401 tips, 1801 internal nodes
\item binindaemonds2007 - 4510 tips, 2108 internal nodes
\item zanne2014 - 31749 tips, 31748 internal nodes
\item smith2011 - 55473 tips, 55338 internal nodes
}
}
\source{
phylocom
}
\description{
Phylogenies to use with phylomatic
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/brranching-package.R
\docType{data}
\name{tpl}
\alias{tpl}
\title{Lookup-table for family, genus, and species names for ThePlantList
gymnosperms}
\format{
A data frame with 23,801 rows and 2 variables:
\itemize{
\item family: family name
\item genus: genus name
}
}
\source{
\url{http://www.theplantlist.org/}
}
\description{
These names are from \url{http://www.theplantlist.org/}, collected
on 2015-11-11, and are from version 1.1 of their data. This data is
used in the function \code{\link[=phylomatic_names]{phylomatic_names()}}
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phylomatic_names.R
\name{phylomatic_names}
\alias{phylomatic_names}
\title{Phylomatic names}
\usage{
phylomatic_names(taxa, format = "isubmit", db = "ncbi", ...)
}
\arguments{
\item{taxa}{quoted tsn number (taxonomic serial number)}

\item{format}{output format, isubmit (you can paste in to the Phylomatic
website), or 'rsubmit' to use in fxn phylomatic_tree}

\item{db}{One of "ncbi", "itis", or "apg". if you use "apg", no HTTP
requests are made (no internet connection needed), whereas if you use
"ncbi" or "itis" you do need an internet connection. IMPORTANT:
see \strong{Authentication} below if using "ncbi".}

\item{...}{curl options passed on to \code{\link[taxize:tax_name]{taxize::tax_name()}}}
}
\value{
string (e.g., "pinaceae/pinus/pinus_contorta"), in Phylomatic
submission format
}
\description{
Get family names to make Phylomatic input object, and
output input string to Phylomatic for use in the function phylomatic
}
\section{Authentication}{

NCBI Entrez doesn't require that you use an API key, but you get
higher rate limit with a key, from 3 to 10 requests per second, so do
get one. Run \code{taxize::use_entrez()} or see
https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/
for instructions.

NCBI API key handling logic is done inside of the \code{taxize} package, used
inside this function.

Save your API key with the name \code{ENTREZ_KEY} as an R option in your
\code{.Rprofile} file, or as environment variables in either your \code{.Renviron}
file or \code{.bash_profile} file, or \code{.zshrc} file (if you use oh-my-zsh) or
similar. See \link{Startup} for help on R options and environment
variables. You cannot pass in your API key in this function.

Remember to restart your R session (and to start a new shell window/tab
if you're using the shell) to take advantage of the new R options
or environment variables.

We strongly recommend using environment variables over R options.

Note that if you don't have an ENTREZ_KEY set, you'll get a message
about it, but only once during each function call. That is, there
can be of these messages per R session, across function calls.
}

\examples{
\dontrun{
mynames <- c("Poa annua", "Salix goodingii", "Helianthus annuus")
phylomatic_names(taxa = mynames, format='rsubmit')
phylomatic_names(mynames, format='rsubmit', db="apg")
phylomatic_names(mynames, format='isubmit', db="ncbi")
phylomatic_names(mynames, format='isubmit', db="apg")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phylomatic_local.R
\name{phylomatic_local}
\alias{phylomatic_local}
\title{Use Phylomatic locally - ideal for large queries}
\usage{
phylomatic_local(
  taxa,
  taxnames = TRUE,
  storedtree = "R20120829",
  db = "apg",
  lowercase = FALSE,
  nodes = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{taxa}{(character) Phylomatic format input of taxa names. required}

\item{taxnames}{If \code{TRUE} (default), we get the family names for you
to attach to your species names to send to Phylomatic API. If \code{FALSE},
you have to provide the strings in the right format. See Details.}

\item{storedtree}{One of R20120829 (Phylomatic tree R20120829 for plants),
smith2011 (Smith 2011, plants), binindaemonds2007 (Bininda-Emonds 2007,
mammals), or zanne2014 (Zanne et al. 2014, plants). Default: R20120829}

\item{db}{(character) One of "ncbi", "itis", or "apg". Default: apg}

\item{lowercase}{(logical) Convert all chars in taxa file to lowercase.
Default: \code{FALSE}}

\item{nodes}{(logical) label all nodes with default names.
Default: \code{FALSE}}

\item{verbose}{(logical) Print messages. Default: \code{TRUE}}
}
\value{
Newick formatted tree as \code{phylo} object
}
\description{
Use Phylomatic locally - ideal for large queries
}
\details{
uses \code{\link[phylocomr:ph_phylomatic]{phylocomr::ph_phylomatic()}} under the hood

This function uses Phylomatic via Phylocom using the
\pkg{phylocomr} package. The interface is slightly different from
\code{\link[=phylomatic]{phylomatic()}}: there's no tree by URL available, and some of the
parameters are not included here.

If you set \code{taxnames = FALSE}, you need to pass in a character
vector, with each element like this example:
\code{"asteraceae/taraxacum/taraxacum_officinale"}, of the form
\code{"family/genus/genus_specfic epithet"}
}
\examples{
\dontrun{
library('ape')

# Input taxonomic names
taxa <- c("Poa annua", "Phlox diffusa", "Helianthus annuus")
(tree <- phylomatic_local(taxa))
plot(collapse.singles(tree), no.margin=TRUE)

taxa <- c("Poa annua", "Collomia grandiflora", "Lilium lankongense",
"Phlox diffusa", "Iteadaphne caudata", "Gagea sarmentosa",
"Helianthus annuus")
(tree <- phylomatic_local(taxa))
plot(collapse.singles(tree), no.margin=TRUE)

library("taxize")
spp <- names_list("species", 500)
length(spp)
(tree <- phylomatic_local(spp))
}
}
