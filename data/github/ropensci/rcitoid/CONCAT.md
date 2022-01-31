rcitoid
=========

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/rcitoid)](https://cranchecks.info/pkgs/rcitoid)
[![R-check](https://github.com/ropensci/rcitoid/workflows/R-check/badge.svg)](https://github.com/ropensci/rcitoid/actions?query=workflow%3AR-check)
[![codecov](https://codecov.io/gh/ropensci/rcitoid/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/rcitoid)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/rcitoid)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/rcitoid)](https://cran.r-project.org/package=rcitoid)




Client for the Citoid service https://www.mediawiki.org/wiki/Citoid

docs: https://en.wikipedia.org/api/rest_v1/#!/Citation/getCitation

There are two functions, both of which do the same things, except:

- `cit_oid()`: parses text
- `cit_oid_()`: does not parse text, you can parse later yourself

Even with `cit_oid()` though, you get a list of lists, and you may
want to parse it to a data.frame. See an example below.

## Install

Stable version


```r
install.packages("rcitoid")
```

Development version


```r
remotes::install_github("ropensci/rcitoid")
```

Load the package



```r
library("rcitoid")
```

## get citation data

use underscore method to get text


```r
cit_oid_("10.1108/jd-12-2013-0166")
#> [[1]]
#> [1] "[{\"key\":\"6FVPATDA\",\"version\":0,\"itemType\":\"journalArticle\",\"tags\":[],\"publicationTitle\":\"Journal of Documentation\",\"journalAbbreviation\":\"Journal of Documentation\",\"volume\":\"71\",\"issue\":\"2\",\"language\":\"en\",\"ISSN\":[\"0022-0418\"],\"date\":\"2015-03-09\",\"pages\":\"253–277\",\"DOI\":\"10.1108/JD-12-2013-0166\",\"url\":\"https://www.emerald.com/insight/content/doi/10.1108/JD-12-2013-0166/full/html\",\"title\":\"Setting our bibliographic references free: towards open citation data\",\"libraryCatalog\":\"DOI.org (Crossref)\",\"accessDate\":\"2020-11-26\",\"shortTitle\":\"Setting our bibliographic references free\",\"author\":[[\"Silvio\",\"Peroni\"],[\"Alexander\",\"Dutton\"],[\"Tanya\",\"Gray\"],[\"David\",\"Shotton\"]],\"source\":[\"Zotero\"]}]"
#> attr(,"type")
#> [1] "json"
```

## get citation data

DOI


```r
cit_oid("10.1108/jd-12-2013-0166")
#> [[1]]
#> [[1]]$key
#> [1] "HCXFCV8F"
#> 
#> [[1]]$version
#> [1] 0
#> 
#> [[1]]$itemType
#> [1] "journalArticle"
#> 
...
```

PMID


```r
cit_oid(30446726)
#> [[1]]
#> [[1]]$key
#> [1] "JDITNMEK"
#> 
#> [[1]]$version
#> [1] 0
#> 
#> [[1]]$itemType
#> [1] "journalArticle"
#> 
...
```

PMCID


```r
cit_oid("PMC4679344")
#> [[1]]
#> [[1]]$key
#> [1] "SY68KLQD"
#> 
#> [[1]]$version
#> [1] 0
#> 
#> [[1]]$itemType
#> [1] "journalArticle"
#> 
...
```

ISBN


```r
cit_oid(1439895619)
#> [[1]]
#> [[1]]$itemType
#> [1] "book"
#> 
#> [[1]]$title
#> [1] "Agroecology : the ecology of sustainable food systems"
#> 
#> [[1]]$oclc
#> [1] "908080219"
#> 
...
```

## parse to data.frame

because the resulting data is nested and can have missing data slots,
it's probably easier to get raw text and manipulate from there.


```r
library(dplyr)

pmid <- c(30446726, 30722046, 30687373, 30688010)
pmcid <- c("PMC4679344", "PMC6347797", "PMC6347793")
isbn <- 1439895619
dois <- c("10.1109/jsac.2011.110806", "10.1007/s00422-006-0078-4",
  "10.5040/9781474219624-0044", "10.1109/icemi.2009.5274826",
  "10.1109/wispnet.2017.8299996")
res <- cit_oid_(id = c(pmid, pmcid, isbn, dois))
tbl_df(bind_rows(lapply(res, jsonlite::fromJSON)))
#> # A tibble: 13 x 33
#>    key   version itemType tags  title pages ISSN  journalAbbrevia…
#>    <chr>   <int> <chr>    <lis> <chr> <chr> <lis> <chr>           
#>  1 P4CD…       0 journal… <df[… Enha… 555–… <chr… Mucosal Immunol 
#>  2 QR7H…       0 journal… <df[… Shar… 1113… <chr… Mol Biol Evol   
#>  3 SLTV…       0 journal… <df[… Resp… 1981  <chr… Front Plant Sci 
#>  4 CWF3…       0 journal… <df[… Mixe… 604–… <chr… Integr Zool     
#>  5 YG5T…       0 journal… <lis… ESMO… 2–30  <chr… Int J Gynecol C…
#>  6 6W95…       0 journal… <lis… Effi… <NA>  <chr… J Orthop Surg R…
#>  7 R3ZP…       0 journal… <lis… Iden… <NA>  <chr… J Hematol Oncol 
#>  8 <NA>       NA book     <NUL… Agro… <NA>  <NUL… <NA>            
#>  9 G24A…       0 journal… <lis… Anti… 1392… <chr… IEEE J. Select.…
#> 10 SXSZ…       0 journal… <lis… The … 193–… <chr… Biol Cybern     
#> 11 LG4G…       0 book     <lis… The … <NA>  <NUL… <NA>            
#> 12 778X…       0 confere… <df[… Desi… 1–47… <NUL… <NA>            
#> 13 P5X4…       0 confere… <lis… Traf… 1412… <NUL… <NA>            
#> # … with 25 more variables: publicationTitle <chr>, volume <chr>, issue <chr>,
#> #   date <chr>, abstractNote <chr>, DOI <chr>, extra <chr>,
#> #   libraryCatalog <chr>, url <chr>, accessDate <chr>, author <list>,
#> #   PMID <chr>, PMCID <chr>, source <list>, shortTitle <chr>, oclc <chr>,
#> #   ISBN <list>, place <chr>, numPages <chr>, contributor <list>,
#> #   language <chr>, publisher <chr>, editor <list>, proceedingsTitle <chr>,
#> #   conferenceName <chr>
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/rcitoid/issues)
* License: MIT
* Get citation information for `rcitoid` in R doing `citation(package = 'rcitoid')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
rcitoid 0.1.0
=============

### NEW FEATURES

* Released to CRAN
## Test environments

* local OS X install, R 3.5.2 patched
* ubuntu 14.04 (on travis-ci), R 3.5.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

* Note about license:
License components with restrictions and base license permitting such:
  MIT + file LICENSE
File 'LICENSE':
  YEAR: 2019
  COPYRIGHT HOLDER: Scott Chamberlain

## Reverse dependencies

This is a new submission, so there are no reverse dependencies.

---

This is a new release. I have read and agree to the the 
CRAN policies at https://cran.r-project.org/web/packages/policies.html

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

Please note that the rcitoid project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.

### See rOpenSci [contributing guide](https://devguide.ropensci.org/contributingguide.html)
for further details.

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.

### Prefer to Email? 

Open an issue instead.
<!-- If authentication is involved: do not share your username/password, or api keys/tokens in this issue - most likely the maintainer will have their own equivalent key -->

<!-- Do not share screen shots of code. Share actual code in text format. -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/). If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
rcitoid
=========

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/rcitoid)](https://cranchecks.info/pkgs/rcitoid)
[![R-check](https://github.com/ropensci/rcitoid/workflows/R-check/badge.svg)](https://github.com/ropensci/rcitoid/actions?query=workflow%3AR-check)
[![codecov](https://codecov.io/gh/ropensci/rcitoid/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/rcitoid)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/rcitoid)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/rcitoid)](https://cran.r-project.org/package=rcitoid)


```{r echo=FALSE}
library("knitr")
hook_output <- knitr::knit_hooks$get("output")
knitr::knit_hooks$set(output = function(x, options) {
   lines <- options$output.lines
   if (is.null(lines)) {
     return(hook_output(x, options))  # pass to default hook
   }
   x <- unlist(strsplit(x, "\n"))
   more <- "..."
   if (length(lines)==1) {        # first n lines
     if (length(x) > lines) {
       # truncate the output, but add ....
       x <- c(head(x, lines), more)
     }
   } else {
     x <- c(if (abs(lines[1])>1) more else NULL,
            x[lines],
            if (length(x)>lines[abs(length(lines))]) more else NULL
           )
   }
   # paste these lines together
   x <- paste(c(x, ""), collapse = "\n")
   hook_output(x, options)
 })

knitr::opts_chunk$set(
  warning = FALSE,
  message = FALSE,
  collapse = TRUE,
  comment = "#>"
)
```

Client for the Citoid service https://www.mediawiki.org/wiki/Citoid

docs: https://en.wikipedia.org/api/rest_v1/#!/Citation/getCitation

There are two functions, both of which do the same things, except:

- `cit_oid()`: parses text
- `cit_oid_()`: does not parse text, you can parse later yourself

Even with `cit_oid()` though, you get a list of lists, and you may
want to parse it to a data.frame. See an example below.

## Install

Stable version

```{r eval=FALSE}
install.packages("rcitoid")
```

Development version

```{r eval=FALSE}
remotes::install_github("ropensci/rcitoid")
```

Load the package


```{r}
library("rcitoid")
```

## get citation data

use underscore method to get text

```{r}
cit_oid_("10.1108/jd-12-2013-0166")
```

## get citation data

DOI

```{r output.lines=1:10}
cit_oid("10.1108/jd-12-2013-0166")
```

PMID

```{r output.lines=1:10}
cit_oid(30446726)
```

PMCID

```{r output.lines=1:10}
cit_oid("PMC4679344")
```

ISBN

```{r output.lines=1:10}
cit_oid(1439895619)
```

## parse to data.frame

because the resulting data is nested and can have missing data slots,
it's probably easier to get raw text and manipulate from there.

```{r}
library(dplyr)

pmid <- c(30446726, 30722046, 30687373, 30688010)
pmcid <- c("PMC4679344", "PMC6347797", "PMC6347793")
isbn <- 1439895619
dois <- c("10.1109/jsac.2011.110806", "10.1007/s00422-006-0078-4",
  "10.5040/9781474219624-0044", "10.1109/icemi.2009.5274826",
  "10.1109/wispnet.2017.8299996")
res <- cit_oid_(id = c(pmid, pmcid, isbn, dois))
tbl_df(bind_rows(lapply(res, jsonlite::fromJSON)))
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/rcitoid/issues)
* License: MIT
* Get citation information for `rcitoid` in R doing `citation(package = 'rcitoid')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rcitoid-package.R
\docType{package}
\name{rcitoid-package}
\alias{rcitoid-package}
\alias{rcitoid}
\title{rcitoid}
\description{
Client for Citoid (\url{https://www.mediawiki.org/wiki/Citoid})
}
\author{
Scott Chamberlain
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cit_oid.R
\name{cit_oid}
\alias{cit_oid}
\alias{cit_oid_}
\title{Get Citoid data}
\usage{
cit_oid(id, format = "mediawiki", accept_language = NULL, ...)

cit_oid_(id, format = "mediawiki", accept_language = NULL, ...)
}
\arguments{
\item{id}{(character) id of an article, DOI, ISBN, PMCID or PMID}

\item{format}{(character) the format to use for the resulting citation data.
one of mediawiki (default), zotero, mediawiki-basefields, bibtex}

\item{accept_language}{(character) for some articles the result depends on
the accept language value, so provide it if localized content is required}

\item{...}{curl options passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\value{
list of lists or character, see
http://opencitations.net/index/coci/api/v1 for explanation of the resulting
columns
}
\description{
Get Citoid data
}
\details{
\code{cit_oid_()} gets raw text (either bibtex or JSON), and \code{cit_oid()}
parses the text as appropriate for the type
}
\examples{
url<-"https://en.wikipedia.org/api/rest_v1/data/citation/mediawiki/30446726"
if (crul::ok(url)) {
  pmid1 <- 30446726
  cit_oid(pmid1)
}

\dontrun{
doi1 <- "10.1108/jd-12-2013-0166"
doi2 <- "10.1371/journal.pone.0058568"
pmid1 <- 30446726
pmcid1 <- "PMC4679344"
isbn1 <- 1439895619

# different formats
cit_oid(doi1)
cit_oid(pmid1, format = "mediawiki")
cit_oid(pmid1, format = "zotero")
cit_oid(pmid1, format = "mediawiki-basefields")
cat(cit_oid(pmid1, format = "bibtex")[[1]])

# PMID example
cit_oid(pmid1, verbose = TRUE)

# ISBN example
cit_oid(isbn1, verbose = TRUE)

# PMCID example
cit_oid(pmcid1)

# set the accept language
x <- cit_oid(pmid1, accept_language = "fr-FR", verbose = TRUE)
x <- cit_oid(doi2, accept_language = "de-DE", verbose = TRUE)

# just get raw text/json
cit_oid_(pmcid1)

# many ids at once
cit_oid(id = c(pmid1, pmcid1, isbn1))
cit_oid_(id = c(pmid1, pmcid1, isbn1))
cit_oid_(id = c(pmid1, pmcid1, isbn1), format = "bibtex")
}
}
\references{
https://en.wikipedia.org/api/rest_v1/#!/Citation/getCitation,
https://www.mediawiki.org/wiki/Citoid
}
