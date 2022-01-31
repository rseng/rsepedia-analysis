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
rdryad
======

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE
)
```

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

```{r eval=FALSE}
install.packages("rdryad")
```

development version:

```{r eval=FALSE}
remotes::install_github("ropensci/rdryad")
```

```{r}
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
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{search_dryad}
\alias{search_dryad}
\title{Search metadata for search terms using regex}
\usage{
search_dryad(...)
}
\description{
This function is defunct
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{dryad_package_dois}
\alias{dryad_package_dois}
\title{Get file DOIs for a Dryad package DOI}
\usage{
dryad_package_dois(doi, ...)
}
\arguments{
\item{...}{ignored}
}
\description{
This function is defunct
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{download_dryadmetadata}
\alias{download_dryadmetadata}
\title{Download metadata for individual Dryad id's}
\usage{
download_dryadmetadata(...)
}
\description{
This function changed name to \code{\link[=dr_get_records]{dr_get_records()}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dryad_download.R
\name{dryad_download}
\alias{dryad_download}
\title{dryad_download}
\usage{
dryad_download(dois, ...)
}
\arguments{
\item{dois}{(character) one or more DOIs, required}

\item{...}{Further args passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\value{
file path for the file
}
\description{
Download datasets by their DOI(s)
}
\examples{
\dontrun{
dryad_download(dois = "10.5061/dryad.f385721n")
dois <- c("10.5061/dryad.f385721n", "10.5061/dryad.7ct1n", "10.5061/dryad.1g626")
dryad_download(dois = dois)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdryad-package.R
\docType{package}
\name{rdryad-package}
\alias{rdryad-package}
\alias{rdryad}
\title{Interface to the Dryad Web services}
\description{
Includes access to Dryad's Solr API, OAI-PMH service, and part of
their REST API.
}
\section{Package API}{


The functions match the three major sets of Dryad API routes for
datasets, fiiles and versions.

Datasets:
\itemize{
\item \code{\link[=dryad_dataset]{dryad_dataset()}}
\item \code{\link[=dryad_datasets]{dryad_datasets()}}
\item \code{\link[=dryad_dataset_versions]{dryad_dataset_versions()}}
}

Files:
\itemize{
\item \code{\link[=dryad_files]{dryad_files()}}
\item \code{\link[=dryad_files_download]{dryad_files_download()}}
}

Versions:
\itemize{
\item \code{\link[=dryad_versions]{dryad_versions()}}
\item \code{\link[=dryad_versions_files]{dryad_versions_files()}}
\item \code{\link[=dryad_versions_download]{dryad_versions_download()}}
}
}

\section{Defunct}{


The Dryad Solr API is no longer being updated, so the functions
that used to work with it are all defunct, see \link{solr-defunct}

The Dryad OAI-PMH service is no longer being updated, so the functions
that used to work with it are all defunct, see \link{oai-defunct}

More defunct functions:
\itemize{
\item \code{\link[=dryad_metadata]{dryad_metadata()}}
\item \code{\link[=dryad_package_dois]{dryad_package_dois()}}
\item \code{\link[=handle2doi]{handle2doi()}}
\item \code{\link[=doi2handle]{doi2handle()}}
\item \code{\link[=dryad_files]{dryad_files()}}
\item \code{\link[=dryad_fetch]{dryad_fetch()}} - use instead \code{\link[=dryad_files_download]{dryad_files_download()}} or
\code{\link[=dryad_versions_download]{dryad_versions_download()}}
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\name{dryad_dataset}
\alias{dryad_dataset}
\title{Get datasets by DOI(s)}
\usage{
dryad_dataset(dois, ...)
}
\arguments{
\item{dois}{(character) one or more DOIs, required}

\item{...}{Further args passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\value{
a list of lists, each named by the input DOI
}
\description{
Get datasets by DOI(s)
}
\examples{
\dontrun{
dryad_dataset(doi = "10.5061/dryad.f385721n")
dois <- c("10.5061/dryad.f385721n", "10.5061/dryad.7ct1n", "10.5061/dryad.1g626")
dryad_dataset(dois = dois)
}
}
\seealso{
Other dryad-datasets: 
\code{\link{dryad_dataset_versions}()},
\code{\link{dryad_datasets}()}
}
\concept{dryad-datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{download_url}
\alias{download_url}
\title{Download url}
\usage{
download_url(...)
}
\description{
This function changed name to \link{dryad_files}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{dryad_metadata}
\alias{dryad_metadata}
\title{Download Dryad file metadata}
\usage{
dryad_metadata(doi, ...)
}
\arguments{
\item{...}{ignored}
}
\description{
This function is defunct
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/files.R
\name{dryad_files}
\alias{dryad_files}
\title{Get metadata information about a file}
\usage{
dryad_files(ids, ...)
}
\arguments{
\item{ids}{(numeric) one or more file ids, required}

\item{...}{Further args passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\value{
a list of lists, each named by the input DOI
}
\description{
Get metadata information about a file
}
\examples{
\dontrun{
dryad_files(ids = 61859)
dryad_files(ids = 61858)
dryad_files(ids = c(61858, 61859))
}
}
\seealso{
Other dryad-files: 
\code{\link{dryad_files_download}()}
}
\concept{dryad-files}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{getalldryad_metadata}
\alias{getalldryad_metadata}
\title{Download metadata for all Dryad oai's for defined time period}
\usage{
getalldryad_metadata(...)
}
\description{
This function is defunct
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{doi2handle}
\alias{doi2handle}
\alias{handle2doi}
\title{Get a Dryad DOI from a handle, and vice versa}
\usage{
doi2handle(...)

handle2doi(...)
}
\arguments{
\item{...}{ignored}
}
\description{
These functions are defunct
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/files.R
\name{dryad_files_download}
\alias{dryad_files_download}
\title{Download a specific file}
\usage{
dryad_files_download(ids, ...)
}
\arguments{
\item{ids}{(numeric) one or more file ids, required}

\item{...}{Further args passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\value{
a list of lists, each named by the input DOI
}
\description{
Download a specific file
}
\note{
UPDATE: we used to not use caching in this fxn; we do now
as of 2020-12-15
}
\examples{
\dontrun{
dryad_files_download(ids = 61858)
dryad_files_download(ids = 61859)
}
}
\seealso{
Other dryad-files: 
\code{\link{dryad_files}()}
}
\concept{dryad-files}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{dryad_fetch}
\alias{dryad_fetch}
\title{Download Dryad files}
\usage{
dryad_fetch(...)
}
\arguments{
\item{...}{ignored}
}
\description{
This function is defunct
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/versions.R
\name{versions}
\alias{versions}
\alias{dryad_versions}
\alias{dryad_versions_files}
\alias{dryad_versions_download}
\title{Get a dataset version by version ID}
\usage{
dryad_versions(ids, ...)

dryad_versions_files(ids, ...)

dryad_versions_download(ids, ...)
}
\arguments{
\item{ids}{(numeric/integer) one or more version ids, required}

\item{...}{Further args passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\value{
a list of lists, each named by the input DOI
}
\description{
Get a dataset version by version ID
}
\details{
\code{dryad_versions()} and \code{dryad_versions_files()}
use async http requests, while \code{dryad_versions_download()}
does not use async
}
\examples{
\dontrun{
dryad_versions(ids = 18774)
dryad_versions_files(ids = 18774)
dryad_versions_download(ids = 18774)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{oai-defunct}
\alias{oai-defunct}
\alias{dr_get_records}
\alias{dr_identify}
\alias{dr_list_records}
\alias{dr_list_identifiers}
\alias{dr_list_metadata_formats}
\alias{dr_list_sets}
\title{Defunct OAI-PMH functions}
\usage{
dr_get_records(...)

dr_identify(...)

dr_list_records(
  prefix = "oai_dc",
  from = NULL,
  until = NULL,
  set = "hdl_10255_3",
  token = NULL,
  as = "df",
  ...
)

dr_list_identifiers(
  prefix = "oai_dc",
  from = NULL,
  until = NULL,
  set = "hdl_10255_3",
  token = NULL,
  as = "df",
  ...
)

dr_list_metadata_formats(...)

dr_list_sets(token = NULL, as = "df", ...)
}
\arguments{
\item{...}{ignored}
}
\description{
Defunct OAI-PMH functions

Download metadata for individual Dryad id's

Learn about the Dryad OAI-PMH service.

List Dryad records

Gets OAI Dryad identifiers

Get available Dryad metadata formats

List the sets in the Dryad metadata repository.
}
\details{
The Dryad OAI-PMH service is no longer being updated
See http://wiki.datadryad.org/Old:Dryad_API#OAI-PMH

Defunct functions:
\itemize{
\item \code{dr_get_records}
\item \code{dr_identify}
\item \code{dr_list_records}
\item \code{dr_list_identifiers}
\item \code{dr_list_metadata_formats}
\item \code{dr_list_sets}
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\name{dryad_dataset_versions}
\alias{dryad_dataset_versions}
\title{Get dataset versions by DOI(s)}
\usage{
dryad_dataset_versions(dois, ...)
}
\arguments{
\item{dois}{(character) one or more DOIs, required}

\item{...}{Further args passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\value{
a list of lists, each named by the input DOI
}
\description{
Get dataset versions by DOI(s)
}
\examples{
\dontrun{
x = dryad_dataset_versions(dois = "10.5061/dryad.f385721n")
x
dois <- c("10.5061/dryad.f385721n", "10.5061/dryad.7ct1n", "10.5061/dryad.1g626")
dryad_dataset_versions(dois = dois)
}
}
\seealso{
Other dryad-datasets: 
\code{\link{dryad_datasets}()},
\code{\link{dryad_dataset}()}
}
\concept{dryad-datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{solr-defunct}
\alias{solr-defunct}
\alias{d_solr_search}
\alias{d_solr_facet}
\alias{d_solr_group}
\alias{d_solr_highlight}
\alias{d_solr_mlt}
\alias{d_solr_stats}
\title{Defunct Solr functions}
\usage{
d_solr_search(...)

d_solr_facet(...)

d_solr_group(...)

d_solr_highlight(...)

d_solr_mlt(...)

d_solr_stats(...)
}
\arguments{
\item{...}{ignored}
}
\description{
Defunct Solr functions

Search the Dryad Solr endpoint
}
\details{
The Dryad Solr service is no longer being updated
See http://wiki.datadryad.org/Old:Dryad_API#SOLR_search_access

Defunct functions:
\itemize{
\item \code{d_solr_search}
\item \code{d_solr_facet}
\item \code{d_solr_group}
\item \code{d_solr_highlight}
\item \code{d_solr_mlt}
\item \code{d_solr_stats}
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\name{dryad_datasets}
\alias{dryad_datasets}
\title{List datasets}
\usage{
dryad_datasets(...)
}
\arguments{
\item{...}{Further args passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\value{
a tibble
}
\description{
List datasets
}
\examples{
\dontrun{
(x <- dryad_datasets())
x$meta
x$links
x$data
}
}
\seealso{
Other dryad-datasets: 
\code{\link{dryad_dataset_versions}()},
\code{\link{dryad_dataset}()}
}
\concept{dryad-datasets}
