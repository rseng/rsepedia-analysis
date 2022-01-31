zbank
=====



[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/zbank)](https://cranchecks.info/pkgs/zbank)
[![R-check](https://github.com/ropensci/zbank/workflows/R-check/badge.svg)](https://github.com/ropensci/zbank/actions?query=workflow%3AR-check)
[![codecov.io](https://codecov.io/github/ropensci/zbank/coverage.svg?branch=master)](https://codecov.io/github/ropensci/zbank?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/zbank)](https://github.com/metacran/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/zbank)](https://cran.r-project.org/package=zbank)

ZooBank API Client

## ZooBank API Docs

See http://zoobank.org/Api

## High vs. Low level package APIs

__High level API__

High level functions do the HTTP request and parse data to a data.frame for
ease of downstream use. The high level functions have no underscore on the end
of the function name, e.g., `zb_name_usages`

__Low level API__

The parsing to data.frame in the high level API does take extra time. The low
level API only does the HTTP request, and gives back JSON without doing any
more parsing. The low level functions DO have an underscore on the end
of the function name, e.g., `zb_name_usages_`

## Install

CRAN version


```r
install.packages("zbank")
```

Development version


```r
remotes::install_github("ropensci/zbank")
```


```r
library("zbank")
```

## Examples

Name usages


```r
zb_name_usages(name = "Pseudanthias carlsoni")
#> # A tibble: 1 x 14
#>   tnuuuid originalreferen… protonymuuid label value lsid  parentname namestring
#>   <chr>   <chr>            <chr>        <chr> <chr> <chr> <chr>      <chr>     
#> 1 6ea8bb… 427d7953-e8fc-4… 6ea8bb2a-a5… carl… carl… urn:… ""         carlsoni  
#> # … with 6 more variables: rankgroup <chr>, usageauthors <chr>,
#> #   taxonnamerankid <chr>, parentusageuuid <chr>, cleanprotonym <chr>,
#> #   nomenclaturalcode <chr>
```

Publications


```r
zb_publications(query = "pyle")
#> # A tibble: 159 x 22
#>    referenceuuid label value authorlist year  title citationdetails volume
#>    <chr>         <chr> <chr> <chr>      <chr> <chr> <chr>           <chr> 
#>  1 a91facc3-f39… <Uns… <Uns… <Unspecif… ""    [Ori… ""              ""    
#>  2 24689ae4-c77… Alon… Alon… Alonso-Za… "201… Manu… "<em>ZooKeys</… "550" 
#>  3 cb272efd-80d… Aran… Aran… Arango, B… "201… Thre… "<em>ZooKeys</… "835" 
#>  4 913bb1fb-1c2… Asgh… Asgh… Asghari, … "201… Desc… "<em>Zootaxa</… "3986"
#>  5 f8ece6ce-e77… Bald… Bald… Baldwin, … "199… <i>B… "<I>Ichthyolog… "45"  
#>  6 72d9641b-2a9… Bass… Bass… Basset, Y… "199… Para… "<em>AMBIO: A … ""    
#>  7 5cb7ee8b-042… Bass… Bass… Basset, Y… "200… Quan… "<em>Bioscienc… "50"  
#>  8 a06a8b51-46b… Bass… Bass… Bassett, … "199… Para… "<em>AMBIO: A … ""    
#>  9 fafe53c6-fdd… Bous… Bous… Boustany,… "200… Sate… "<em>Nature</e… "415" 
#> 10 18985db2-cc5… Bush… Bush… Bush, Eli… "199… Taxo… "<em>Harvard P… "2"   
#> # … with 149 more rows, and 14 more variables: number <chr>, edition <chr>,
#> #   publisher <chr>, placepublished <chr>, pagination <chr>, startpage <chr>,
#> #   endpage <chr>, language <chr>, languageid <chr>, referencetype <chr>,
#> #   lsid <chr>, parentreferenceid <chr>, parentreference <chr>, authors <list>
```

Authors


```r
zb_authors(query = "Schmutz")
#> # A tibble: 3 x 9
#>   agentnameid label value zblsid familyname givenname preferreduuid agentid
#>   <chr>       <chr> <chr> <chr>  <chr>      <chr>     <chr>         <chr>  
#> 1 F16D374C-5… Achi… Achi… "urn:… Achitte-S… Helga C.  F16D374C-531… F16D37…
#> 2 EC923CC6-E… Schm… Schm… "urn:… Schmutz    Karl      EC923CC6-E5E… EC923C…
#> 3 2E04A84F-4… Schm… Schm… ""     Schmutzler Clarence… 2E04A84F-459… 2E04A8…
#> # … with 1 more variable: isuser <chr>
```

Get info by any ZooBank identifier


```r
zb_id(id = "6EA8BB2A-A57B-47C1-953E-042D8CD8E0E2")
#> # A tibble: 3 x 10
#>   identifier identifierdomain abbreviation identifierurl registeringagen…
#>   <chr>      <chr>            <chr>        <chr>         <chr>           
#> 1 66491      CAS Ichthy Spec… CAS_SPC      http://resea… ""              
#> 2 643345     Taxonomic Seria… ITIS-TSN     http://www.i… ""              
#> 3 urn:lsid:… ZooBank Nomencl… ZB-Act       http://zooba… "William N."    
#> # … with 5 more variables: registeringagentfamilyname <chr>,
#> #   registeringagentorganizationname <chr>, identifieruuid <chr>,
#> #   domainlogourl <chr>, resolutionnote <chr>
```

Matching taxon name service


```r
zb_matching(id = "FFF7160A-372D-40E9-9611-23AF5D9EAC4C")
#> # A tibble: 36 x 8
#>    protonymuuid acceptedprotony… protonymuuidarr… fullnamestring
#>    <chr>        <chr>            <chr>            <chr>         
#>  1 FFF7160A-37… FFF7160A-372D-4… FBDC898C-F1EA-4… Cheilodipteru…
#>  2 FFF7160A-37… FFF7160A-372D-4… 97CE20CD-4D6E-4… Gasterosteus …
#>  3 FFF7160A-37… FFF7160A-372D-4… 97CE20CD-4D6E-4… Gaſteroſteus …
#>  4 FFF7160A-37… FFF7160A-372D-4… 4D4EC609-D241-4… Pomatomus sal…
#>  5 FFF7160A-37… FFF7160A-372D-4… 4D4EC609-D241-4… Pomatomus sal…
#>  6 FFF7160A-37… FFF7160A-372D-4… 67D18558-763F-4… Temnodon salt…
#>  7 FFF7160A-37… FFF7160A-372D-4… FBDC898C-F1EA-4… Cheilodipteru…
#>  8 FFF7160A-37… FFF7160A-372D-4… 97CE20CD-4D6E-4… Gasterosteus …
#>  9 FFF7160A-37… FFF7160A-372D-4… 97CE20CD-4D6E-4… Gaſteroſteus …
#> 10 FFF7160A-37… FFF7160A-372D-4… 4D4EC609-D241-4… Pomatomus sal…
#> # … with 26 more rows, and 4 more variables: nomenclaturalcodeid <int>,
#> #   taxonrank <chr>, synonymtype <int>, referencecount <int>
```

ZooBank usage stats


```r
zb_stats(start_date = "2018-03-01", end_date = "2018-04-01")
#> # A tibble: 96 x 3
#>    identifierdomain day        recordcount
#>    <chr>            <chr>      <chr>      
#>  1 ZooBank Author   2018-03-01 39         
#>  2 ZooBank Author   2018-03-02 25         
#>  3 ZooBank Author   2018-03-03 9          
#>  4 ZooBank Author   2018-03-04 6          
#>  5 ZooBank Author   2018-03-05 23         
#>  6 ZooBank Author   2018-03-06 17         
#>  7 ZooBank Author   2018-03-07 26         
#>  8 ZooBank Author   2018-03-08 20         
#>  9 ZooBank Author   2018-03-09 19         
#> 10 ZooBank Author   2018-03-10 3          
#> # … with 86 more rows
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/zbank/issues).
* License: MIT
* Get citation information for `zbank` in R doing `citation(package = 'zbank')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![rofooter](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
zbank 0.1.0
===========

### NEW FEATURES

* released to CRAN
## Test environments

* local OS X install, R 3.5.1 patched
* ubuntu 14.04 (on travis-ci), R 3.5.1
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

  License components with restrictions and base license permitting such:
    MIT + file LICENSE
  File 'LICENSE':
    YEAR: 2018
    COPYRIGHT HOLDER: Scott Chamberlain

## Reverse dependencies

This is a new submission, so there are no reverse dependencies.

--------

This is a new release. I have read and agree to the the CRAN policies at https://cran.r-project.org/web/packages/policies.html

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

* Submit an issue on the [Issues page](https://github.com/ropensci/zbank/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/zbank.git`
* Make sure to track progress upstream (i.e., on our version of `zbank` at `ropensci/zbank`) by doing `git remote add upstream https://github.com/ropensci/zbank.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Push up to your account
* Submit a pull request to home base at `ropensci/zbank`

### Also, check out our [discussion forum](https://discuss.ropensci.org)
<!-- If authentication is involved: do not share your username/password, or api keys/tokens in this issue - most likely the maintainer will have their own equivalent key -->

<!-- Do not share screen shots of code. Share actual code in text format. -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/). If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
zbank
=====

```{r echo=FALSE}
library("knitr")
library("zbank")
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

[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/zbank)](https://cranchecks.info/pkgs/zbank)
[![R-check](https://github.com/ropensci/zbank/workflows/R-check/badge.svg)](https://github.com/ropensci/zbank/actions?query=workflow%3AR-check)
[![codecov.io](https://codecov.io/github/ropensci/zbank/coverage.svg?branch=master)](https://codecov.io/github/ropensci/zbank?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/zbank)](https://github.com/metacran/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/zbank)](https://cran.r-project.org/package=zbank)

ZooBank API Client

## ZooBank API Docs

See http://zoobank.org/Api

## High vs. Low level package APIs

__High level API__

High level functions do the HTTP request and parse data to a data.frame for
ease of downstream use. The high level functions have no underscore on the end
of the function name, e.g., `zb_name_usages`

__Low level API__

The parsing to data.frame in the high level API does take extra time. The low
level API only does the HTTP request, and gives back JSON without doing any
more parsing. The low level functions DO have an underscore on the end
of the function name, e.g., `zb_name_usages_`

## Install

CRAN version

```{r eval=FALSE}
install.packages("zbank")
```

Development version

```{r eval=FALSE}
remotes::install_github("ropensci/zbank")
```

```{r}
library("zbank")
```

## Examples

Name usages

```{r}
zb_name_usages(name = "Pseudanthias carlsoni")
```

Publications

```{r}
zb_publications(query = "pyle")
```

Authors

```{r}
zb_authors(query = "Schmutz")
```

Get info by any ZooBank identifier

```{r}
zb_id(id = "6EA8BB2A-A57B-47C1-953E-042D8CD8E0E2")
```

Matching taxon name service

```{r}
zb_matching(id = "FFF7160A-372D-40E9-9611-23AF5D9EAC4C")
```

ZooBank usage stats

```{r}
zb_stats(start_date = "2018-03-01", end_date = "2018-04-01")
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/zbank/issues).
* License: MIT
* Get citation information for `zbank` in R doing `citation(package = 'zbank')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![rofooter](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zb_publications.R
\name{zb_publications}
\alias{zb_publications}
\alias{zb_publications_}
\title{Publications}
\usage{
zb_publications(id = NULL, query = NULL, parse = TRUE, ...)

zb_publications_(id, query, ...)
}
\arguments{
\item{id}{(integer/numeric) A publication identifier}

\item{query}{(character) Query terms}

\item{parse}{(logical) Whether to parse to list (\code{FALSE}) or
data.frame (\code{TRUE}). If \code{TRUE}, we also give back a
tibble for easier consumption. Default: \code{TRUE}}

\item{...}{Curl options passed to \code{\link[crul]{HttpClient}}}
}
\description{
Publications
}
\examples{
\dontrun{
zb_publications(id = "427D7953-E8FC-41E8-BEA7-8AE644E6DE77")
zb_publications(query = "pyle")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zb_stats.R
\name{zb_stats}
\alias{zb_stats}
\alias{zb_stats_}
\title{Get statistics on Zoobank activity}
\usage{
zb_stats(start_date, end_date, period = "day", parse = TRUE, ...)

zb_stats_(start_date, end_date, period, ...)
}
\arguments{
\item{start_date}{(date/character) a start date}

\item{end_date}{(date/character) an end date}

\item{period}{(character) the period. Default: day}

\item{parse}{(logical) Whether to parse to list (\code{FALSE}) or
data.frame (\code{TRUE}). If \code{TRUE}, we also give back a
tibble for easier consumption. Default: \code{TRUE}}

\item{...}{Curl options passed to \code{\link[crul]{HttpClient}}}
}
\description{
Get statistics on Zoobank activity
}
\examples{
\dontrun{
zb_stats(start_date = "2018-03-01", end_date = "2018-04-01")
zb_stats(start_date = "2018-03-01", end_date = "2018-04-01", 
 parse = FALSE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zb_authors.R
\name{zb_authors}
\alias{zb_authors}
\alias{zb_authors_}
\title{Search for authors}
\usage{
zb_authors(id = NULL, query = NULL, parse = TRUE, ...)

zb_authors_(id = NULL, query = NULL, ...)
}
\arguments{
\item{id}{(integer/numeric) An author identifier}

\item{query}{(character) Query terms}

\item{parse}{(logical) Whether to parse to list (\code{FALSE}) or
data.frame (\code{TRUE}). If \code{TRUE}, we also give back a
tibble for easier consumption. Default: \code{TRUE}}

\item{...}{Curl options passed to \code{\link[crul]{HttpClient}}}
}
\description{
Search for authors
}
\examples{
\dontrun{
zb_authors(id = "8C466CBE-3F7D-4DC9-8CBD-26DD3F57E212")
zb_authors(query = "Schmutz")
zb_authors(query = "Pyle")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zb_matching.R
\name{zb_matching}
\alias{zb_matching}
\title{Matching taxon name service}
\usage{
zb_matching(id, parse = TRUE, ...)
}
\arguments{
\item{id}{(integer/numeric) Any ZooBank identifier, for taxon, author
or publication. required}

\item{parse}{(logical) Whether to parse to list (\code{FALSE}) or
data.frame (\code{TRUE}). If \code{TRUE}, we also give back a
tibble for easier consumption. Default: \code{TRUE}}

\item{...}{Curl options passed to \code{\link[crul]{HttpClient}}}
}
\description{
Matching taxon name service
}
\examples{
\dontrun{
zb_matching(id = "FFF7160A-372D-40E9-9611-23AF5D9EAC4C")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zb_name_usages.R
\name{zb_name_usages}
\alias{zb_name_usages}
\alias{zb_name_usages_}
\title{Name usages}
\usage{
zb_name_usages(name = NULL, id = NULL, query = NULL, parse = TRUE, ...)

zb_name_usages_(name = NULL, id = NULL, query = NULL, ...)
}
\arguments{
\item{name}{(character) A taxonomic name}

\item{id}{(integer/numeric) A taxonomic identifier}

\item{query}{(character) A taxonomic name to query}

\item{parse}{(logical) Whether to parse to list (\code{FALSE}) or
data.frame (\code{TRUE}). If \code{TRUE}, we also give back a
tibble for easier consumption. Default: \code{TRUE}}

\item{...}{Curl options passed to \code{\link[crul]{HttpClient}}}
}
\description{
Name usages
}
\examples{
\dontrun{
zb_name_usages(name = "Pseudanthias carlsoni")
zb_name_usages(id = "6EA8BB2A-A57B-47C1-953E-042D8CD8E0E2")
zb_name_usages(query = "pyle")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zb_id.R, R/zb_matching.R
\name{zb_id}
\alias{zb_id}
\alias{zb_id_}
\alias{zb_matching_}
\title{Get data by identifier}
\usage{
zb_id(id, parse = TRUE, ...)

zb_id_(id, ...)

zb_matching_(id, ...)
}
\arguments{
\item{id}{(integer/numeric) Any ZooBank identifier, for taxon, author
or publication. required}

\item{parse}{(logical) Whether to parse to list (\code{FALSE}) or
data.frame (\code{TRUE}). If \code{TRUE}, we also give back a
tibble for easier consumption. Default: \code{TRUE}}

\item{...}{Curl options passed to \code{\link[crul]{HttpClient}}}
}
\description{
Get data by identifier
}
\examples{
\dontrun{
zb_id(id = "6EA8BB2A-A57B-47C1-953E-042D8CD8E0E2")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zbank-package.R
\docType{package}
\name{zbank-package}
\alias{zbank-package}
\alias{zbank}
\title{zbank}
\description{
ZooBank Client
}
\section{ZooBank API Docs}{

See http://zoobank.org/Api
}

\author{
Scott Chamberlain
}
\keyword{package}
