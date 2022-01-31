rdatacite
=========



[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/rdatacite)](https://cranchecks.info/pkgs/rdatacite)
[![R-CMD-check](https://github.com/ropensci/rdatacite/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/rdatacite/actions?query=workflow%3AR-CMD-check)
[![codecov.io](https://codecov.io/github/ropensci/rdatacite/coverage.svg?branch=master)](https://codecov.io/github/ropensci/rdatacite?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/rdatacite)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/rdatacite)](https://cran.r-project.org/package=rdatacite)
[![DOI](https://zenodo.org/badge/2521192.svg)](https://zenodo.org/badge/latestdoi/2521192)

`rdatacite` provides programmatic accesses to DataCite (https://datacite.org/) metadata

* REST API. Docs: https://support.datacite.org/docs/api and https://support.datacite.org/reference

`rdatacite` docs: https://docs.ropensci.org/rdatacite

Package API:

 - `dc_providers`
 - `dc_reports`
 - `dc_check`
 - `dc_events`
 - `dc_dois`
 - `dc_clients`
 - `dc_client_prefixes`
 - `dc_provider_prefixes`
 - `dc_status`
 - `dc_prefixes`
 - `dc_activities`

## Installation

Stable CRAN version


```r
install.packages("rdatacite")
```

Development version from github


```r
pak::pkg_install("ropensci/rdatacite")
```


```r
library('rdatacite')
```

## Result objects

Outputs from nearly all `rdatacite` functions will be of class `dc`, an S3 class that's 
simply a named list of results. You can easily remove the class via `unclass()`.
The `print.dc` method prints the data.frame for the `data`, `included`, and `reports`
slots if they exist, but hides the `meta` named list. You can get to the metadata by
indexing to it like `$meta`.

## Searching

You may want to start with `dc_dois()`.


```r
dc_dois(query = "climate change")
#> datacite: dois
#> found: 85075, pages: 400, page: 1
#> slots: data, meta, links
#> $data
#> # A tibble: 25 x 4
#>    id    type  attributes$doi $identifiers $creators $titles $publisher
#>    <chr> <chr> <chr>          <list>       <list>    <list>  <chr>     
#>  1 10.1… dois  10.15786/20.5… <list [0]>   <df[,6] … <df[,1… Mountain …
#>  2 10.2… dois  10.25675/1021… <list [0]>   <df[,3] … <df[,2… Mountain …
#>  3 10.2… dois  10.25675/1021… <list [0]>   <df[,3] … <df[,2… Mountain …
#>  4 10.2… dois  10.25675/1021… <list [0]>   <df[,3] … <df[,2… Mountain …
#>  5 10.2… dois  10.25675/1021… <list [0]>   <df[,3] … <df[,2… Mountain …
#>  6 10.2… dois  10.25676/1112… <list [0]>   <df[,6] … <df[,1… Mountain …
#>  7 10.2… dois  10.25676/1112… <list [0]>   <df[,6] … <df[,1… Mountain …
#>  8 10.2… dois  10.25675/1021… <list [0]>   <df[,6] … <df[,2… Mountain …
#>  9 10.2… dois  10.25675/1021… <list [0]>   <df[,6] … <df[,1… Mountain …
#> 10 10.2… dois  10.25675/1021… <list [0]>   <df[,6] … <df[,1… Mountain …
#> # … with 15 more rows, and 42 more variables: $container <df[,0]>,
#> #   $publicationYear <int>, $subjects <list>, $contributors <list>,
#> #   $dates <list>, $language <chr>, $types$ris <chr>, $$bibtex <chr>,
#> #   $$citeproc <chr>, $$schemaOrg <chr>, $$resourceType <chr>,
#> #   $$resourceTypeGeneral <chr>, $relatedIdentifiers <list>, $sizes <list>,
#> #   $formats <list>, $version <lgl>, $rightsList <list>, $descriptions <list>,
#> #   $geoLocations <list>, $fundingReferences <list>, $url <chr>,
#> #   $contentUrl <lgl>, $metadataVersion <int>, $schemaVersion <chr>,
#> #   $source <chr>, $isActive <lgl>, $state <chr>, $reason <lgl>,
#> #   $viewCount <int>, $downloadCount <int>, $referenceCount <int>,
#> #   $citationCount <int>, $partCount <int>, $partOfCount <int>,
#> #   $versionCount <int>, $versionOfCount <int>, $created <chr>,
#> #   $registered <chr>, $published <lgl>, $updated <chr>,
#> #   relationships$client$data$id <chr>, $$$type <chr>
#> 
#> $included
#> NULL
```

The `query` parameter supports Elasticearch query string queries. Some examples:


```r
# search within a field
dc_dois(query = "publicationYear:2016")
# fuzzy search (via *) on a nested field
dc_dois(query = "creators.familyName:mil*")
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/rdatacite/issues).
* License: MIT
* Get citation information for `rdatacite` in R doing `citation(package = 'rdatacite')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
rdatacite 0.5.2
===============

### MINOR IMPROVEMENTS

* fix breaking test on one of the cran checks (#30)


rdatacite 0.5.0
===============

### NEW FEATURES

* Major refactor to work with the new DataCite API: all functions from the previous version are defunct; all OAI-PMH functions are gone; new functions all start with `dc_` (#24) (#29)

### MINOR IMPROVEMENTS

* all examples check if DataCite API is up before running (#28)


rdatacite 0.4.2
===============

### MINOR IMPROVEMENTS

* fix to two fixtures that had non-ascii text in them, that were causing tests to fail (#25)


rdatacite 0.4.0
===============

### MINOR IMPROVEMENTS

* pagination fixes (#18)
* fix unused httr package warning, flagged by cran team (#21)
* add .github PR and issue templates


rdatacite 0.3.0
===============

### NEW FEATURES

* Gains new functions for working with the DataCite REST API:
`dc_data_center`, `dc_data_centers`, `dc_member`, `dc_members`,
`dc_work`, `dc_works` (#13)
* Now using new version of solrium package - users shouldn't see any differences (#16)

### BUG FIXES

* Fix scientific notation (#15)
* Fix `vapply` error (#14)



rdatacite 0.1.0
===============

### NEW FEATURES

* Released to CRAN.
## Test environments

* local OS X install, R 3.6.3 RC
* ubuntu 14.04 (on travis-ci), R 3.6.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are no reverse dependencies.

---

This version fixes a problem in the test suite causing a R CMD CHECK failure on R devel linux.

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

* Submit an issue on the [Issues page](https://github.com/ropensci/rdatacite/issues)

### Issues and Pull Requests

If you are considering a pull request, you may want to open an issue first to discuss with the maintainer(s).

### Code contributions

* Fork this repo to your GitHub account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/rdatacite.git`
* Make sure to track progress upstream (i.e., on our version of `rdatacite` at `ropensci/rdatacite`) by doing `git remote add upstream https://github.com/ropensci/rdatacite.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch - see <https://guides.github.com/introduction/flow/> for how to contribute by branching, making changes, then submitting a pull request)
* Push up to your account
* Submit a pull request to home base (likely master branch, but check to make sure) at `ropensci/rdatacite`

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
rdatacite
=========

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE
)
```

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/rdatacite)](https://cranchecks.info/pkgs/rdatacite)
[![R-CMD-check](https://github.com/ropensci/rdatacite/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/rdatacite/actions?query=workflow%3AR-CMD-check)
[![codecov.io](https://codecov.io/github/ropensci/rdatacite/coverage.svg?branch=master)](https://codecov.io/github/ropensci/rdatacite?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/rdatacite)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/rdatacite)](https://cran.r-project.org/package=rdatacite)
[![DOI](https://zenodo.org/badge/2521192.svg)](https://zenodo.org/badge/latestdoi/2521192)

`rdatacite` provides programmatic accesses to DataCite (https://datacite.org/) metadata

* REST API. Docs: https://support.datacite.org/docs/api and https://support.datacite.org/reference

`rdatacite` docs: https://docs.ropensci.org/rdatacite

Package API:

```{r echo=FALSE, comment=NA, results='asis'}
cat(paste(" -", paste(sprintf("`%s`", getNamespaceExports("rdatacite")), collapse = "\n - ")))
```

## Installation

Stable CRAN version

```{r eval=FALSE}
install.packages("rdatacite")
```

Development version from github

```{r eval=FALSE}
pak::pkg_install("ropensci/rdatacite")
```

```{r}
library('rdatacite')
```

## Result objects

Outputs from nearly all `rdatacite` functions will be of class `dc`, an S3 class that's 
simply a named list of results. You can easily remove the class via `unclass()`.
The `print.dc` method prints the data.frame for the `data`, `included`, and `reports`
slots if they exist, but hides the `meta` named list. You can get to the metadata by
indexing to it like `$meta`.

## Searching

You may want to start with `dc_dois()`.

```{r}
dc_dois(query = "climate change")
```

The `query` parameter supports Elasticearch query string queries. Some examples:

```{r eval=FALSE}
# search within a field
dc_dois(query = "publicationYear:2016")
# fuzzy search (via *) on a nested field
dc_dois(query = "creators.familyName:mil*")
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/rdatacite/issues).
* License: MIT
* Get citation information for `rdatacite` in R doing `citation(package = 'rdatacite')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dc_providers.R
\name{dc_providers}
\alias{dc_providers}
\title{DataCite REST API: providers}
\usage{
dc_providers(
  ids = NULL,
  query = NULL,
  year = NULL,
  region = NULL,
  organization_type = NULL,
  focus_area = NULL,
  include = NULL,
  limit = 25,
  page = 1,
  cursor = NULL,
  ...
)
}
\arguments{
\item{ids}{(character) one or more provider IDs}

\item{query}{(character) query string}

\item{year}{(character) year}

\item{region}{(character) region name}

\item{organization_type}{(character) organization type}

\item{focus_area}{(character) focus area}

\item{include}{(character) vector of fields to return}

\item{limit}{(numeric/integer) results per page}

\item{page}{(numeric/integer) result page, the record to start at}

\item{cursor}{(character) page cursor (used instead of \code{limit} param)}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
DataCite REST API: providers
}
\examples{
\dontrun{
if (dc_check()) {
x <- dc_providers()
x
dc_providers(limit = 3)
dc_providers(ids = x$data$id[1:5])
}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dc_status.R
\name{dc_check}
\alias{dc_check}
\title{check if the DataCite API is up or not}
\usage{
dc_check(...)
}
\value{
boolean
}
\description{
check if the DataCite API is up or not
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dc_provider_prefixes.R
\name{dc_provider_prefixes}
\alias{dc_provider_prefixes}
\title{DataCite REST API: provider prefixes}
\usage{
dc_provider_prefixes(include = NULL, limit = 25, page = 1, cursor = NULL, ...)
}
\arguments{
\item{include}{(character) vector of fields to return}

\item{limit}{(numeric/integer) results per page}

\item{page}{(numeric/integer) result page, the record to start at}

\item{cursor}{(character) page cursor (used instead of \code{limit} param)}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
DataCite REST API: provider prefixes
}
\examples{
\dontrun{
if (dc_check()) {
x <- dc_provider_prefixes()
x
dc_provider_prefixes(limit = 3)
}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dc_clients.R
\name{dc_clients}
\alias{dc_clients}
\title{DataCite REST API: clients}
\usage{
dc_clients(
  ids = NULL,
  query = NULL,
  year = NULL,
  provider_id = NULL,
  software = NULL,
  include = NULL,
  limit = 25,
  page = 1,
  cursor = NULL,
  ...
)
}
\arguments{
\item{ids}{(character) one or more client IDs}

\item{query}{(character) Query string}

\item{year}{(integer/numeric/character) a year}

\item{provider_id}{a provider ID}

\item{software}{no idea what should go here, anyone?}

\item{include}{(character) vector of fields to return}

\item{limit}{(numeric/integer) results per page}

\item{page}{(numeric/integer) the page to get results for. default: 1}

\item{cursor}{(character) page cursor (used instead of \code{limit} param).
to use cursor pagination, set \code{cursor = 1}, then use the link in
\verb{$links$next}}

\item{...}{curl options passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\description{
DataCite REST API: clients
}
\examples{
\dontrun{
if (dc_check()) {
x <- dc_clients()
x
dc_clients(x$data$id[1])
dc_clients(x$data$id[1:2], verbose = TRUE)
}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dc_cn.R
\name{dc_cn}
\alias{dc_cn}
\title{DataCite content negotation}
\usage{
dc_cn(dois, format = "bibtex", style = "apa", locale = "en-US", ...)
}
\arguments{
\item{dois}{(character) one or more DOIs}

\item{format}{Name of the format. One of "rdf-xml", "turtle",
"citeproc-json", "schemaorg", "codemeta", "text", "ris", "bibtex"
(default), "datacite-xml", "datacite-json", "bibentry", or
"jats".}

\item{style}{a CSL style (for text format only). See
‘rcrossref::get_styles()’ for options. Default: 'apa'. If there's
a style that DataCite doesn't support you'll get a
(500) Internal Server Error}

\item{locale}{Language locale. See ‘?Sys.getlocale’}

\item{...}{curl options passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\description{
DataCite content negotation
}
\examples{
\dontrun{
dc_cn("10.5281/zenodo.50213")
dc_cn(c("10.5281/zenodo.50213", "10.5281/zenodo.57081"), "text")
dc_cn(c("a-bad-doi", "10.5281/zenodo.50213", "10.5281/zenodo.57081"), "text")
}
}
\references{
https://support.datacite.org/docs/datacite-content-resolver
}
\seealso{
see also \code{rcrossref::cr_cn} for a more general purpose
content negotation interface
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dc_status.R
\name{dc_status}
\alias{dc_status}
\title{DataCite REST API: status of the API}
\usage{
dc_status(...)
}
\arguments{
\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
DataCite REST API: status of the API
}
\examples{
\dontrun{
if (dc_check()) {
dc_status()
}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dc_reports.R
\name{dc_reports}
\alias{dc_reports}
\title{DataCite REST API: reports}
\usage{
dc_reports(
  ids = NULL,
  platform = NULL,
  report_name = NULL,
  report_id = NULL,
  release = NULL,
  created = NULL,
  created_by = NULL,
  include = NULL,
  limit = 25,
  page = 1,
  ...
)
}
\arguments{
\item{ids}{(character) one or more report IDs}

\item{platform}{(character) Name of the Platform the usage is being
requested for. This can be omitted if the service provides usage for
only one platform.}

\item{report_name}{(character) The long name of the report}

\item{report_id}{(character) The report ID or code or shortname. Typically
this will be the same code provided in the Report parameter of the request}

\item{release}{(character) The release or version of the report}

\item{created}{(character) Time the report was prepared. Format as defined
by date-time - RFC3339}

\item{created_by}{(character) Name of the organization producing the report}

\item{include}{(character) vector of fields to return}

\item{limit}{(numeric/integer) results per page}

\item{page}{(numeric/integer) result page, the record to start at}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
DataCite REST API: reports
}
\examples{
\dontrun{
if (dc_check()) {
x <- dc_reports()
x
dc_reports(created = "2019-08-01T07:00:00.000Z")
dc_reports(created_by = "urn:node:GOA")
dc_reports(limit = 3)
# dc_reports(ids = x$reports$id[1:3]) # FIXME: doesn't work
}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dc_events.R
\name{dc_events}
\alias{dc_events}
\title{DataCite REST API: events}
\usage{
dc_events(
  ids = NULL,
  query = NULL,
  subj_id = NULL,
  obj_id = NULL,
  doi = NULL,
  orcid = NULL,
  prefix = NULL,
  subtype = NULL,
  subject = NULL,
  source_id = NULL,
  registrant_id = NULL,
  relation_type_id = NULL,
  issn = NULL,
  publication_year = NULL,
  year_month = NULL,
  include = NULL,
  sort = NULL,
  limit = 25,
  page = 1,
  cursor = NULL,
  ...
)
}
\arguments{
\item{ids}{(character) one or more event IDs}

\item{query}{(character) Query for any event information}

\item{subj_id}{(character) The identifier for the event subject, expressed
as a URL. For example: \verb{https://doi.org/10.7272/q6qn64nk}}

\item{obj_id}{(character) The identifier for the event object, expressed
as a URL. For example: \verb{https://doi.org/10.7272/q6qn64nk}}

\item{doi}{(character) The subj-id or obj-id of the event, expressed as
a DOI. For example: \code{10.7272/q6qn64nk}}

\item{orcid}{(character) an ORCID, presumably}

\item{prefix}{(character) The DOI prefix of the subj-id or obj-id of the
event. For example: \code{10.7272}}

\item{subtype}{(character) xxx}

\item{subject}{(character) xxx}

\item{source_id}{(character) a source ID. See Details}

\item{registrant_id}{(character)}

\item{relation_type_id}{(character) a relation-type ID. See Details}

\item{issn}{(character) an ISSN, presumably}

\item{publication_year}{(character) the publication year}

\item{year_month}{(character) The year and month in which the event
occurred, in the format \code{YYYY-MM}. For example \code{2018-08}}

\item{include}{(character) vector of fields to return}

\item{sort}{(character) variable to sort by}

\item{limit}{(numeric/integer) results per page}

\item{page}{(numeric/integer) the page to get results for. default: 1}

\item{cursor}{(character) page cursor (used instead of \code{limit} param).
to use cursor pagination, set \code{cursor = 1}, then use the link in
\verb{$links$next}}

\item{...}{curl options passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\description{
DataCite REST API: events
}
\details{
See https://support.datacite.org/docs/eventdata-guide for
details on possible values for parameters
}
\examples{
\dontrun{
if (dc_check()) {
# dc_events(query = "birds")
}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dc_dois.R
\name{dc_dois}
\alias{dc_dois}
\title{DataCite REST API: dois}
\usage{
dc_dois(
  ids = NULL,
  query = NULL,
  created = NULL,
  registered = NULL,
  provider_id = NULL,
  client_id = NULL,
  person_id = NULL,
  resource_type_id = NULL,
  subject = NULL,
  schema_version = NULL,
  random = NULL,
  sample_size = NULL,
  sample_group = NULL,
  include = NULL,
  sort = NULL,
  limit = 25,
  page = 1,
  cursor = NULL,
  ...
)
}
\arguments{
\item{ids}{(character) one or more DOIs}

\item{query}{(character) Query string. See Querying below.}

\item{created}{(character) metadata where year of DOI creation is \code{created}.
See Filtering Responses below.}

\item{registered}{(character) metadata where year of DOI registration
is \code{year}. See Filtering Responses below.}

\item{provider_id}{(character) metadata associated with a specific DataCite
provider. See Filtering Responses below.}

\item{client_id}{(character) metadata associated with a specific DataCite
client. See Filtering Responses below.}

\item{person_id}{(character) metadata associated with a specific person's
ORCID iD. See Filtering Responses below.}

\item{resource_type_id}{(character) metadata for a specific
resourceTypeGeneral. See Filtering Responses below.}

\item{subject}{(character)}

\item{schema_version}{(character) metadata where schema version of the
deposited metadata is \code{schema-version}. See Filtering Responses below.}

\item{random}{(logical) return random set of results, can be combined
with any kind of query. default: \code{FALSE}.}

\item{sample_size}{(character)}

\item{sample_group}{(character)}

\item{include}{(character) vector of fields to return}

\item{sort}{(character) variable to sort by}

\item{limit}{(numeric/integer) results per page}

\item{page}{(numeric/integer) the page to get results for. default: 1}

\item{cursor}{(character) page cursor (used instead of \code{limit} param).
to use cursor pagination, set \code{cursor = 1}, then use the link in
\verb{$links$next}}

\item{...}{curl options passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\description{
DataCite REST API: dois
}
\section{Querying}{

See https://support.datacite.org/docs/api-queries for details
}

\section{Filtering Responses}{

See
https://support.datacite.org/docs/api-queries#section-filtering-list-responses
for details
}

\examples{
\dontrun{
if (dc_check()) {
x <- dc_dois()
x
dc_dois(query = "birds")
dc_dois(query = "climate change")
dc_dois(query = "publicationYear:2016")
x <- dc_dois(query = "creators.familyName:mil*", verbose = TRUE)
lapply(x$data$attributes$creators, "[[", "familyName")
x <- dc_dois(query = "titles.title:climate +change")
lapply(x$data$attributes$titles, "[[", "title")
dc_dois(client_id = "dryad.dryad")
dc_dois(x$data$id[1])
dc_dois(x$data$id[1:3])
dc_dois("10.5281/zenodo.1308060")

# pagination
dc_dois(limit = 1)
x <- dc_dois(cursor = 1)
x$links$`next`
}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dc_activities.R
\name{dc_activities}
\alias{dc_activities}
\title{DataCite REST API: activities}
\usage{
dc_activities(
  ids = NULL,
  query = NULL,
  limit = 25,
  page = 1,
  cursor = NULL,
  ...
)
}
\arguments{
\item{ids}{(character) one or more activity IDs}

\item{query}{(character) Query string}

\item{limit}{(numeric/integer) results per page}

\item{page}{(numeric/integer) the page to get results for. default: 1}

\item{cursor}{(character) page cursor (used instead of \code{limit} param).
to use cursor pagination, set \code{cursor = 1}, then use the link in
\verb{$links$next}}

\item{...}{curl options passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\description{
DataCite REST API: activities
}
\details{
for more info on the \verb{/activities} route see
https://support.datacite.org/docs/tracking-provenance
}
\examples{
\dontrun{
if (dc_check()) {
x <- dc_activities()
x
# dc_activities(x$data$id[1]) # FIXME: doesn't work, returns no data
# dc_activities(query = "ecology") # FIXME: this thlimit a 500 error
}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdatacite-package.R
\docType{package}
\name{rdatacite-package}
\alias{rdatacite-package}
\alias{rdatacite}
\title{rdatacite}
\description{
DataCite R client
}
\section{HTTP Requests}{

All HTTP requests are GET requests, and are sent with the following
headers:
\itemize{
\item \verb{Accept: application/vnd.api+json; version=2}
\item \verb{User-Agent: r-curl/4.3 crul/0.9.0 rOpenSci(rdatacite/0.5.0)}
\item \verb{X-USER-AGENT: r-curl/4.3 crul/0.9.0 rOpenSci(rdatacite/0.5.0)}
}

The user-agent strings change as the versions of each package change.
}

\section{Methods in the package}{

\itemize{
\item \code{\link[=dc_providers]{dc_providers()}}
\item \code{\link[=dc_reports]{dc_reports()}}
\item \code{\link[=dc_check]{dc_check()}}
\item \code{\link[=dc_events]{dc_events()}}
\item \code{\link[=dc_dois]{dc_dois()}}
\item \code{\link[=dc_clients]{dc_clients()}}
\item \code{\link[=dc_client_prefixes]{dc_client_prefixes()}}
\item \code{\link[=dc_provider_prefixes]{dc_provider_prefixes()}}
\item \code{\link[=dc_status]{dc_status()}}
\item \code{\link[=dc_prefixes]{dc_prefixes()}}
\item \code{\link[=dc_activities]{dc_activities()}}
}
}

\section{rdatacite defunct functions}{

\itemize{
\item \code{dc_data_center}
\item \code{dc_data_centers}
\item \code{dc_facet}
\item \code{dc_member}
\item \code{dc_members}
\item \code{dc_mlt}
\item \code{dc_oai_getrecord}
\item \code{dc_oai_identify}
\item \code{dc_oai_listidentifiers}
\item \code{dc_oai_listmetadataformats}
\item \code{dc_oai_listrecords}
\item \code{dc_oai_listsets}
\item \code{dc_search}
\item \code{dc_stats}
\item \code{dc_work}
\item \code{dc_works}
}
}

\section{Content negotation}{

For content negotation see \code{rcrossref::cr_cn()}, which can be used for
Crossref, DataCite and Medra DOIs
}

\section{GraphGL API}{

rdatacite does not support the GraphGL API
https://support.datacite.org/docs/datacite-graphql-api-guide - we suggest
trying the \code{ghql} package (https://github.com/ropensci/ghql/)
}

\author{
Scott Chamberlain \email{myrmecocystus@gmail.com}
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dc_client_prefixes.R
\name{dc_client_prefixes}
\alias{dc_client_prefixes}
\title{DataCite REST API: client prefixes}
\usage{
dc_client_prefixes(
  query = NULL,
  year = NULL,
  client_id = NULL,
  prefix_id = NULL,
  sort = NULL,
  include = NULL,
  limit = 25,
  page = 1,
  cursor = NULL,
  ...
)
}
\arguments{
\item{query}{(character) Query string}

\item{year}{(integer/numeric/character) a year}

\item{client_id}{a client ID}

\item{prefix_id}{a prefix ID}

\item{sort}{(character) variable to sort by}

\item{include}{(character) vector of fields to return}

\item{limit}{(numeric/integer) results per page}

\item{page}{(numeric/integer) the page to get results for. default: 1}

\item{cursor}{(character) page cursor (used instead of \code{limit} param).
to use cursor pagination, set \code{cursor = 1}, then use the link in
\verb{$links$next}}

\item{...}{curl options passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\description{
DataCite REST API: client prefixes
}
\examples{
\dontrun{
if (dc_check()) {
x <- dc_client_prefixes()
x
}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dc_prefixes.R
\name{dc_prefixes}
\alias{dc_prefixes}
\title{DataCite REST API: prefixes}
\usage{
dc_prefixes(include = NULL, limit = 25, page = 1, cursor = NULL, ...)
}
\arguments{
\item{include}{(character) vector of fields to return}

\item{limit}{(numeric/integer) results per page}

\item{page}{(numeric/integer) result page, the record to start at}

\item{cursor}{(character) page cursor (used instead of \code{limit} param)}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
DataCite REST API: prefixes
}
\examples{
\dontrun{
if (dc_check()) {
x <- dc_prefixes()
x
dc_prefixes(limit = 3)
}}
}
