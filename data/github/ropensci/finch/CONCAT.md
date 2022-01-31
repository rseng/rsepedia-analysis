finch
=====



[![R-check](https://github.com/ropensci/finch/workflows/R-check/badge.svg)](https://github.com/ropensci/finch/actions?query=workflow%3AR-check)
[![cran checks](https://cranchecks.info/badges/worst/finch)](https://cranchecks.info/pkgs/finch)
[![codecov](https://codecov.io/gh/ropensci/finch/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/finch)
[![cran version](https://www.r-pkg.org/badges/version/finch)](https://cran.r-project.org/package=finch)

`finch` parses Darwin Core simple and archive files

Docs: <https://docs.ropensci.org/finch/>

* Darwin Core description at Biodiversity Information Standards site <http://rs.tdwg.org/dwc.htm>
* Darwin Core at Wikipedia <https://en.wikipedia.org/wiki/Darwin_Core>

## Install

Stable version


```r
install.packages("finch")
```

Development version, from GitHub


```r
remotes::install_github("ropensci/finch")
```


```r
library("finch")
```

## Parse

To parse a simple darwin core file like

```
<?xml version="1.0" encoding="UTF-8"?>
<SimpleDarwinRecordSet
 xmlns="http://rs.tdwg.org/dwc/xsd/simpledarwincore/"
 xmlns:dc="http://purl.org/dc/terms/"
 xmlns:dwc="http://rs.tdwg.org/dwc/terms/"
 xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
 xsi:schemaLocation="http://rs.tdwg.org/dwc/xsd/simpledarwincore/ ../../xsd/tdwg_dwc_simple.xsd">
 <SimpleDarwinRecord>
  <dwc:occurrenceID>urn:catalog:YPM:VP.057488</dwc:occurrenceID>
  <dc:type>PhysicalObject</dc:type>
  <dc:modified>2009-02-12T12:43:31</dc:modified>
  <dc:language>en</dc:language>
  <dwc:basisOfRecord>FossilSpecimen</dwc:basisOfRecord>
  <dwc:institutionCode>YPM</dwc:institutionCode>
  <dwc:collectionCode>VP</dwc:collectionCode>
  <dwc:catalogNumber>VP.057488</dwc:catalogNumber>
  <dwc:individualCount>1</dwc:individualCount>
  <dwc:locationID xsi:nil="true"/>
  <dwc:continent>North America</dwc:continent>
  <dwc:country>United States</dwc:country>
  <dwc:countryCode>US</dwc:countryCode>
  <dwc:stateProvince>Montana</dwc:stateProvince>
  <dwc:county>Garfield</dwc:county>
  <dwc:scientificName>Tyrannosourus rex</dwc:scientificName>
  <dwc:genus>Tyrannosourus</dwc:genus>
  <dwc:specificEpithet>rex</dwc:specificEpithet>
  <dwc:earliestPeriodOrHighestSystem>Creataceous</dwc:earliestPeriodOrHighestSystem>
  <dwc:latestPeriodOrHighestSystem>Creataceous</dwc:latestPeriodOrHighestSystem>
  <dwc:earliestEonOrHighestEonothem>Late Cretaceous</dwc:earliestEonOrHighestEonothem>
  <dwc:latestEonOrHighestEonothem>Late Cretaceous</dwc:latestEonOrHighestEonothem>
 </SimpleDarwinRecord>
</SimpleDarwinRecordSet>
```

This file is in this package as an example file, get the file, then `simple()`


```r
file <- system.file("examples", "example_simple_fossil.xml", package = "finch")
out <- simple_read(file)
```

Index to `meta`, `dc` or `dwc`


```r
out$dc
#> [[1]]
#> [[1]]$type
#> [1] "PhysicalObject"
#> 
#> 
#> [[2]]
#> [[2]]$modified
#> [1] "2009-02-12T12:43:31"
#> 
#> 
#> [[3]]
#> [[3]]$language
#> [1] "en"
```

## Parse Darwin Core Archive

To parse a Darwin Core Archive like can be gotten from GBIF use `dwca_read()`

There's an example Darwin Core Archive:


```r
file <- system.file("examples", "0000154-150116162929234.zip", package = "finch")
(out <- dwca_read(file, read = TRUE))
#> <gbif dwca>
#>   Package ID: 6cfaaf9c-d518-4ca3-8dc5-f5aadddc0390
#>   No. data sources: 10
#>   No. datasets: 3
#>   Dataset occurrence.txt: [225 X 443]
#>   Dataset multimedia.txt: [15 X 1]
#>   Dataset verbatim.txt: [209 X 443]
```

List files in the archive


```r
out$files
#> $xml_files
#> [1] "/Users/sckott/Library/Caches/R/finch/0000154-150116162929234/meta.xml"    
#> [2] "/Users/sckott/Library/Caches/R/finch/0000154-150116162929234/metadata.xml"
#> 
#> $txt_files
#> [1] "/Users/sckott/Library/Caches/R/finch/0000154-150116162929234/citations.txt" 
#> [2] "/Users/sckott/Library/Caches/R/finch/0000154-150116162929234/multimedia.txt"
#> [3] "/Users/sckott/Library/Caches/R/finch/0000154-150116162929234/occurrence.txt"
#> [4] "/Users/sckott/Library/Caches/R/finch/0000154-150116162929234/rights.txt"    
#> [5] "/Users/sckott/Library/Caches/R/finch/0000154-150116162929234/verbatim.txt"  
...
```

High level metadata for the whole archive


```r
out$emlmeta
#> additionalMetadata:
#>   metadata:
#>     gbif:
#>       citation:
#>         identifier: 0000154-150116162929234
#>         citation: GBIF Occurrence Download 0000154-150116162929234
#>       physical:
#>         objectName: []
#>         characterEncoding: UTF-8
#>         dataFormat:
#>           externallyDefinedFormat:
#>             formatName: Darwin Core Archive
#>         distribution:
#>           online:
#>             url:
#>               function: download
#>               url: http://api.gbif.org/v1/occurrence/download/request/0000154-150116162929234.zip
#> dataset:
#>   title: GBIF Occurrence Download 0000154-150116162929234
#>   creator:
...
```

High level metadata for each data file, there's many files, but we'll just look at one


```r
hm <- out$highmeta
head( hm$occurrence.txt )
#>   index                                        term delimitedBy
#> 1     0         http://rs.gbif.org/terms/1.0/gbifID        <NA>
#> 2     1           http://purl.org/dc/terms/abstract        <NA>
#> 3     2       http://purl.org/dc/terms/accessRights        <NA>
#> 4     3      http://purl.org/dc/terms/accrualMethod        <NA>
#> 5     4 http://purl.org/dc/terms/accrualPeriodicity        <NA>
#> 6     5      http://purl.org/dc/terms/accrualPolicy        <NA>
```

You can get the same metadata as above for each dataset that went into the tabular dataset downloaded


```r
out$dataset_meta[[1]]
```

View one of the datasets, brief overview.


```r
head( out$data[[1]][,c(1:5)] )
#>      gbifID abstract accessRights accrualMethod accrualPeriodicity
#> 1  50280003       NA                         NA                 NA
#> 2 477550574       NA                         NA                 NA
#> 3 239703844       NA                         NA                 NA
#> 4 239703843       NA                         NA                 NA
#> 5 239703833       NA                         NA                 NA
#> 6 477550692       NA                         NA                 NA
```

You can also give `dwca()` a local directory, or url that contains a Darwin Core Archive.

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/finch/issues).
* License: MIT
* Get citation information for `finch` in R doing `citation(package = 'finch')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![rofooter](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
finch 0.4.0
===========

### BUG FIXES

* fix various package imports: `plyr` and `rappdirs` no longer needed; import `EML::read_eml` instead of the whole package; import `digest::digest` instead of whole package; import `hoardr::hoard` instead of whole package (#27)


finch 0.3.0
===========

### BUG FIXES

* fix to unit tests and `dwc_read()` for a new version of `EML` package (v0.2) (#26) thanks @cboettig


finch 0.2.0
===========

### CACHING CHANGES

Caching has changed in `finch`. We changed to using package `hoardr` for managing caching. Now with `hoardr` on the package loading we create an object that holds methods and info about where to cache, with operating specific routes. In addition, you can set your own cache directory (and we do this in examples/tests using a temp dir instead of user dir). 

The old functions `dwca_cache_delete`, `dwca_cache_delete_all`, `dwca_cache_details`, and `dwca_cache_list` are defunct and replaced with the single `dwca_cache` object. The `dwca_cache` object is an `R6` object that has methods/functions and variables. See the `?dwca_cache` and `?finch-defunct` manual files for details.

### BUG FIXES

* fix to `dwca_read()`: `...` wasn't being passed on to `data.table::fread` internally (#18) (#19) thanks @gustavobio !

### MINOR IMPROVEMENTS

* replaced Suggested package `httr` with `crul` (#16)
* using markdown docs now (#24)
* improvement to `dwca_read()` to download zip files as binary because without that wasn't working on Windows machines (#17) thanks @gustavobio !
* fix for failures on CRAN: don't write to user directory in examples/tests (#23)


finch 0.1.0
===========

### NEW FEATURES

* Released to CRAN.
## Test environments

* local OS X install, R 4.0.2 patched
* ubuntu 14.04 (on travis-ci), R 4.0.2
* win-builder (devel and release)

## R CMD check results

No problems found

## Reverse dependencies

There are 2 reverse dependencies - no problems were found.

---

This version fixes package import problems.

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

### Bugs?

* Submit an issue on the [Issues page](https://github.com/ropensci/finch/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/finch.git`
* Make sure to track progress upstream (i.e., on our version of `finch` at `ropensci/finch`) by doing `git remote add upstream https://github.com/ropensci/finch.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Please do write a test(s) for your changes if they affect code and not just docs
* Push up to your account
* Submit a pull request to home base at `ropensci/finch`

### Check out our [discussion forum](https://discuss.ropensci.org)
<!-- If authentication is involved: do not share your username/password, or api keys/tokens in this issue - most likely the maintainer will have their own equivalent key -->

<!-- Do not share screen shots of code. Share actual code in text format. -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/). If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
# Platform

|field    |value                                       |
|:--------|:-------------------------------------------|
|version  |R version 4.0.2 Patched (2020-06-30 r78761) |
|os       |macOS Catalina 10.15.6                      |
|system   |x86_64, darwin17.0                          |
|ui       |X11                                         |
|language |(EN)                                        |
|collate  |en_US.UTF-8                                 |
|ctype    |en_US.UTF-8                                 |
|tz       |US/Pacific                                  |
|date     |2020-08-10                                  |

# Dependencies

|package |old   |new   |Δ  |
|:-------|:-----|:-----|:--|
|finch   |0.3.0 |0.4.0 |*  |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*finch
=====

```{r echo=FALSE}
hook_output <- knitr::knit_hooks$get("output")
knitr::knit_hooks$set(output = function(x, options) {
  lines <- options$output.lines
  if (is.null(lines)) {
    return(hook_output(x, options))  # pass to default hook
  }
  x <- unlist(strsplit(x, "\n"))
  more <- "..."
  if (length(lines) == 1) {        # first n lines
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

[![R-check](https://github.com/ropensci/finch/workflows/R-check/badge.svg)](https://github.com/ropensci/finch/actions?query=workflow%3AR-check)
[![cran checks](https://cranchecks.info/badges/worst/finch)](https://cranchecks.info/pkgs/finch)
[![codecov](https://codecov.io/gh/ropensci/finch/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/finch)
[![cran version](https://www.r-pkg.org/badges/version/finch)](https://cran.r-project.org/package=finch)

`finch` parses Darwin Core simple and archive files

Docs: <https://docs.ropensci.org/finch/>

* Darwin Core description at Biodiversity Information Standards site <http://rs.tdwg.org/dwc.htm>
* Darwin Core at Wikipedia <https://en.wikipedia.org/wiki/Darwin_Core>

## Install

Stable version

```{r eval=FALSE}
install.packages("finch")
```

Development version, from GitHub

```{r eval=FALSE}
remotes::install_github("ropensci/finch")
```

```{r}
library("finch")
```

## Parse

To parse a simple darwin core file like

```
<?xml version="1.0" encoding="UTF-8"?>
<SimpleDarwinRecordSet
 xmlns="http://rs.tdwg.org/dwc/xsd/simpledarwincore/"
 xmlns:dc="http://purl.org/dc/terms/"
 xmlns:dwc="http://rs.tdwg.org/dwc/terms/"
 xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
 xsi:schemaLocation="http://rs.tdwg.org/dwc/xsd/simpledarwincore/ ../../xsd/tdwg_dwc_simple.xsd">
 <SimpleDarwinRecord>
  <dwc:occurrenceID>urn:catalog:YPM:VP.057488</dwc:occurrenceID>
  <dc:type>PhysicalObject</dc:type>
  <dc:modified>2009-02-12T12:43:31</dc:modified>
  <dc:language>en</dc:language>
  <dwc:basisOfRecord>FossilSpecimen</dwc:basisOfRecord>
  <dwc:institutionCode>YPM</dwc:institutionCode>
  <dwc:collectionCode>VP</dwc:collectionCode>
  <dwc:catalogNumber>VP.057488</dwc:catalogNumber>
  <dwc:individualCount>1</dwc:individualCount>
  <dwc:locationID xsi:nil="true"/>
  <dwc:continent>North America</dwc:continent>
  <dwc:country>United States</dwc:country>
  <dwc:countryCode>US</dwc:countryCode>
  <dwc:stateProvince>Montana</dwc:stateProvince>
  <dwc:county>Garfield</dwc:county>
  <dwc:scientificName>Tyrannosourus rex</dwc:scientificName>
  <dwc:genus>Tyrannosourus</dwc:genus>
  <dwc:specificEpithet>rex</dwc:specificEpithet>
  <dwc:earliestPeriodOrHighestSystem>Creataceous</dwc:earliestPeriodOrHighestSystem>
  <dwc:latestPeriodOrHighestSystem>Creataceous</dwc:latestPeriodOrHighestSystem>
  <dwc:earliestEonOrHighestEonothem>Late Cretaceous</dwc:earliestEonOrHighestEonothem>
  <dwc:latestEonOrHighestEonothem>Late Cretaceous</dwc:latestEonOrHighestEonothem>
 </SimpleDarwinRecord>
</SimpleDarwinRecordSet>
```

This file is in this package as an example file, get the file, then `simple()`

```{r}
file <- system.file("examples", "example_simple_fossil.xml", package = "finch")
out <- simple_read(file)
```

Index to `meta`, `dc` or `dwc`

```{r}
out$dc
```

## Parse Darwin Core Archive

To parse a Darwin Core Archive like can be gotten from GBIF use `dwca_read()`

There's an example Darwin Core Archive:

```{r}
file <- system.file("examples", "0000154-150116162929234.zip", package = "finch")
(out <- dwca_read(file, read = TRUE))
```

List files in the archive

```{r output.lines=1:10}
out$files
```

High level metadata for the whole archive

```{r output.lines=1:20}
out$emlmeta
```

High level metadata for each data file, there's many files, but we'll just look at one

```{r}
hm <- out$highmeta
head( hm$occurrence.txt )
```

You can get the same metadata as above for each dataset that went into the tabular dataset downloaded

```{r eval=FALSE}
out$dataset_meta[[1]]
```

View one of the datasets, brief overview.

```{r}
head( out$data[[1]][,c(1:5)] )
```

You can also give `dwca()` a local directory, or url that contains a Darwin Core Archive.

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/finch/issues).
* License: MIT
* Get citation information for `finch` in R doing `citation(package = 'finch')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![rofooter](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simple.R
\name{simple_read}
\alias{simple_read}
\title{Parse a DarwinRecordSet and SimpleDarwinRecordSet files}
\usage{
simple_read(file)
}
\arguments{
\item{file}{(character) A path to a single simple Darwin Core
file in XML format. Required.}
}
\value{
a S3 class \code{dwc_recordset} when a DarwinRecordSet is given, or
a \code{dwc_simplerecordset} when a SimpleDarwinRecordSet is given. In
each case the object is really just a list, with lightweight S3 class
attached for easy downstream usage. Prints summary to screen by default
}
\description{
Parse a DarwinRecordSet and SimpleDarwinRecordSet files
}
\details{
Make sure when reading a DarwinRecordSet to access the chunks by
position rather than name since duplicate names are allowed in chunks.
}
\examples{
\dontrun{
# SimpleDarwinRecordSet examples
file <- system.file("examples", "example_simple.xml", package = "finch")
simple_read(file)
file <- system.file("examples", "example_simple_fossil.xml",
  package = "finch")
simple_read(file)

# DarwinRecordSet examples
file <- system.file("examples", "example_classes_observation.xml",
  package = "finch")
simple_read(file)

file <- system.file("examples", "example_classes_specimen.xml",
  package = "finch")
simple_read(file)

# access elements of the object
file <- system.file("examples", "example_classes_specimen.xml",
  package = "finch")
res <- simple_read(file)
## namespaces
res$meta
## locations
res$locations
## chunks, the first one
res$chunks[[1]]
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{finch-defunct}
\alias{finch-defunct}
\title{Defunct functions in finch}
\description{
\itemize{
\item \code{dwca_cache_delete}: Defunt - see \link{dwca_cache}
\item \code{dwca_cache_delete_all}: Defunt - see \link{dwca_cache}
\item \code{dwca_cache_details}: Defunt - see \link{dwca_cache}
\item \code{dwca_cache_list}: Defunt - see \link{dwca_cache}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/finch-package.R
\docType{package}
\name{finch-package}
\alias{finch-package}
\alias{finch}
\title{finch}
\description{
Parse Darwin Core Archive files
}
\author{
Scott Chamberlain \email{myrmecocystus@gmail.com}
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dwca.R
\name{dwca_read}
\alias{dwca_read}
\title{Parse Darwin Core Archive}
\usage{
dwca_read(input, read = FALSE, ...)
}
\arguments{
\item{input}{(character) Path to local zip file, directory, or a url.
If a URL it must be for a zip file.}

\item{read}{(logical) Whether or not to read in data files. If \code{FALSE},
we give back paths to files only. Default: \code{FALSE}}

\item{...}{Further args passed on to \code{\link[data.table:fread]{data.table::fread()}}}
}
\description{
Parse Darwin Core Archive
}
\details{
Note that sometimes file reads fail. We use \code{\link[data.table:fread]{data.table::fread()}}
internally, which is very fast, but can fail sometimes. If so, try reading
in the data manually.

When you pass in a URL, we use \pkg{rappdirs} to determine cache path, and
if you pass the same URL again, and your cache is not cleared, we'll
pull from the cache. Passing a file or directory on your local system
won't invoke the caching route, but will go directly to the file/directory.
}
\examples{
\dontrun{
# set up a temporary directory for the example
dwca_cache$cache_path_set(path = "finch", type = "tempdir")

dir <- system.file("examples", "0000154-150116162929234", package = "finch")

# Don't read data in
(x <- dwca_read(dir, read=FALSE))
x$files
x$highmeta
x$dataset_meta[[1]]
x$data

# Read data
(x <- dwca_read(dir, read=TRUE))
head(x$data[[1]])

# Can pass in a zip file
zip <- system.file("examples", "0000154-150116162929234.zip",
  package = "finch")
(out <- dwca_read(zip))
out$files
out$highmeta
out$emlmeta
out$dataset_meta

# Can pass in zip file as a url
url <-
"https://github.com/ropensci/finch/blob/master/inst/examples/0000154-150116162929234.zip?raw=true"
(out <- dwca_read(url))

# another url
url <- "http://ipt.jbrj.gov.br/jbrj/archive.do?r=redlist_2013_taxons&v=3.12"
(out <- dwca_read(url))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cache.R
\name{dwca_cache}
\alias{dwca_cache}
\title{Caching}
\description{
Manage cached \code{finch} files with package \pkg{hoardr}
}
\details{
The dafault cache directory is
\code{paste0(rappdirs::user_cache_dir(), "/R/finch")}, but you can set
your own path using \code{cache_path_set()}

\code{cache_delete} only accepts one file name, while
\code{cache_delete_all} doesn't accept any names, but deletes all files.
For deleting many specific files, use \code{cache_delete} in a \code{\link[=lapply]{lapply()}}
type call
}
\section{Useful user functions}{

\itemize{
\item \code{dwca_cache$cache_path_get()} get cache path
\item \code{dwca_cache$cache_path_set()} set cache path
\item \code{dwca_cache$list()} returns a character vector of full
path file names
\item \code{dwca_cache$files()} returns file objects with metadata
\item \code{dwca_cache$details()} returns files with details
\item \code{dwca_cache$delete()} delete specific files
\item \code{dwca_cache$delete_all()} delete all files, returns nothing
}
}

\examples{
\dontrun{
dwca_cache

# list files in cache
dwca_cache$list()

# delete certain database files
# dwca_cache$delete("file path")
# dwca_cache$list()

# delete all files in cache
# dwca_cache$delete_all()
# dwca_cache$list()

# set a different cache path from the default
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{dwca_cache_delete}
\alias{dwca_cache_delete}
\title{This function is defunct.}
\usage{
dwca_cache_delete(...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{dwca_cache_details}
\alias{dwca_cache_details}
\title{This function is defunct.}
\usage{
dwca_cache_details(...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{dwca_cache_delete_all}
\alias{dwca_cache_delete_all}
\title{This function is defunct.}
\usage{
dwca_cache_delete_all(...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{dwca_cache_list}
\alias{dwca_cache_list}
\title{This function is defunct.}
\usage{
dwca_cache_list(...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as.location.R
\name{as.location}
\alias{as.location}
\alias{as.location.character}
\alias{as.location.location}
\alias{print.location}
\title{Convert a path or URL to a location object}
\usage{
as.location(x, ...)

\method{as.location}{character}(x, ...)

\method{as.location}{location}(x, ...)

\method{print}{location}(x, ...)
}
\arguments{
\item{x}{Input, a path or URL}

\item{...}{Ignored.}
}
\description{
Convert a path or URL to a location object
}
\examples{
# A zip file
file <- system.file("examples/0000154-150116162929234.zip",
  package = "finch")
as.location(file)

# A directory
dir <- system.file("examples/0000154-150116162929234",
  package = "finch")
as.location(dir)

# A URL
as.location("https://httpbin.org/get")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dwca_validate.R
\name{dwca_validate}
\alias{dwca_validate}
\title{Validate a Darwin Core Archive}
\usage{
dwca_validate(x, ifModifiedSince = NULL, browse = FALSE, ...)
}
\arguments{
\item{x}{(character) A url for a Darwin Core Archive. If you have a local
Darwin Core Archive, put it up online somewhere. Required.}

\item{ifModifiedSince}{(character) An optional ISO date (yyyy-mm-dd) to
enable conditional get requests, validating archives only if they have
been modified since the given date. This feature requires the archive url
to honor the if-modified-since http header. Apache webservers for example
do this out of the box for static files, but if you use dynamic scripts
to generate the archive on the fly this might not be recognised. Optional.}

\item{browse}{(logical) Browse to generated report or not.
Default: \code{FALSE}}

\item{...}{Curl options passed to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Validate a Darwin Core Archive
}
\details{
Uses the GBIF DCA validator (http://tools.gbif.org/dwca-validator/)
}
\examples{
\dontrun{
x <- "http://rs.gbif.org/datasets/german_sl.zip"
dwca_validate(x)
}
}
