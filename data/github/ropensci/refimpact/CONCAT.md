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
(http:contributor-covenant.org), version 1.0.0, available at 
http://contributor-covenant.org/version/1/0/0/

<!-- README.md is generated from README.Rmd. Please edit that file -->
refimpact
=========

[![Build Status](https://travis-ci.org/ropensci/refimpact.svg?branch=master)](https://travis-ci.org/ropensci/refimpact) 
[![Build status](https://ci.appveyor.com/api/projects/status/jxj1yela4a6ym6wb/branch/master?svg=true)](https://ci.appveyor.com/project/perrystephenson/refimpact/branch/master) 
[![codecov](https://codecov.io/gh/ropensci/refimpact/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/refimpact) 
[![](https://badges.ropensci.org/78_status.svg)](https://github.com/ropensci/onboarding/issues/78) 
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/refimpact)](https://CRAN.R-project.org/package=refimpact)

**refimpact** provides an API wrapper for the UK Research Excellence Framework 2014 Impact Case Studies Database. You can find more information about this database at <http://impact.ref.ac.uk/CaseStudies/>.

The data may be of interest to you if you are interested in:

-   text mining
-   directed graphs
-   policies for research funding

Case studies in the database are licenced under a CC-BY 4.0 license. The full license can be found [here](https://creativecommons.org/licenses/by/4.0/legalcode) and a more user-friendly version of the license can be be obtained [here](https://creativecommons.org/licenses/by/4.0/).

Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.

Installation
------------

### Install from CRAN

``` r
install.packages("refimpact")
```

### Install from Github

``` r
install.packages("devtools")
devtools::install_github("perrystephenson/refimpact")
```

Usage
-----

See the vignette:

``` r
vignette("refimpact")
```

More Information
----------------

For more information about a specific function you can use the help commands (for example `?ref_get`).

To raise bug reports and issues, please use the issue tracker in Github.

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# refimpact 1.0.0

* Major breaking changes. There is now a single user-facing function `ref_get()`
  which takes the API method as an argument. This standardises a lot of the 
  input validation and error handling, as well as reducing the risk of bugs (as 
  there are less lines of code). Functions from previous version of the package
  are still available, but deprecated.
* A vignette has been added, and the help documentation has been improved.
* The package now uses **httr** when calling the API, which improves reliability
  and provides much better error messages to the end-user when things go wrong.
* The `phrase` parameter to the SearchCaseStudies method now allows text queries
  of any length and complexity
* Bundled a `ref_tags` dataset with the package, to save the end-user from 
  having to iterate through the ListTagValues API method in order to find tags
  to use as parameters when searching the database
* Added a contributor code of conduct
* The entire package was re-written.

# refimpact 0.1.0

* Initial release.

## Test environments

* local OS X install, R 3.4.1
* Ubuntu 12.04 (on travis-ci), R 3.4.0, R 3.3.3, R-devel.
* Windows (on AppVeyor), R 3.4.1, R 3.3.3, R-devel

## R CMD check results

0 ERRORs | 0 WARNINGs | 0 NOTEs

devtools::win_builder() shows three potential spelling errors - each of these 
has been checked and confirmed accurate.

## Downstream dependencies

There are currently no downstream dependencies for this package.
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

# refimpact

[![Build Status](https://travis-ci.org/ropensci/refimpact.svg?branch=master)](https://travis-ci.org/ropensci/refimpact)
[![Build status](https://ci.appveyor.com/api/projects/status/jxj1yela4a6ym6wb/branch/master?svg=true)](https://ci.appveyor.com/project/perrystephenson/refimpact/branch/master)
[![codecov](https://codecov.io/gh/ropensci/refimpact/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/refimpact)
[![](https://badges.ropensci.org/78_status.svg)](https://github.com/ropensci/onboarding/issues/78)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/refimpact)](https://CRAN.R-project.org/package=refimpact)



**refimpact** provides an API wrapper for the UK Research Excellence Framework
2014 Impact Case Studies Database. You can find more information about this
database at
[http://impact.ref.ac.uk/CaseStudies/](http://impact.ref.ac.uk/CaseStudies/).

The data may be of interest to you if you are interested in:

- text mining
- directed graphs
- policies for research funding

Case studies in the database are licenced under a CC-BY 4.0 license. The full
license can be found 
[here](https://creativecommons.org/licenses/by/4.0/legalcode) and a more 
user-friendly version of the license can be be obtained 
[here](https://creativecommons.org/licenses/by/4.0/).

Please note that this project is released with a 
[Contributor Code of Conduct](CONDUCT.md). By participating in this project you 
agree to abide by its terms.

## Installation

### Install from CRAN

```{r, eval=FALSE}
install.packages("refimpact")
```

### Install from Github

```{r, eval=FALSE}
install.packages("devtools")
devtools::install_github("perrystephenson/refimpact")
```

## Usage

See the vignette:

```{r, eval=FALSE}
vignette("refimpact")
```

## More Information

For more information about a specific function you can use the help commands 
(for example `?ref_get`). 

To raise bug reports and issues, please use the issue tracker in Github.

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)

---
title: "UK REF Impact Case Studies"
author: "Perry Stephenson"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{refimpact}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r message=FALSE}
library(refimpact)
```

## Introduction

This package is an API wrapper around the REF Impact Case Studies database API. 
Chances are that if you're looking at this package, you already know what this
dataset is, and you probably know roughly what you're looking for. 

If you have stumbled upon this package however, and you want to know more about
the dataset, you can head [here](http://impact.ref.ac.uk) to find out more. If
you are thinking of using this dataset as a toy dataset for learning, then you 
might find this dataset useful for text mining, amongst other things.

## Core functions

The core function for this package is `ref_get()`, which takes an API method as
the first argument, and some optional arguments depending on the method.

The API methods available are detailed below, but presented here for quick 
reference:

* SearchCaseStudies
* ListUnitsOfAssessment
* ListTagTypes
* ListTagValues
* ListInstitutions

## SearchCaseStudies

This is the core method of the API, and the most important for users of this
package. The search method requires a compulsory argument to the `ref_get()`
function: `query`. This argument takes a list of query parameters, which can be
as simple as a single Case Study ID, which returns a single record. A query
returning a single record is shown below to demonstrate the syntax and the 
returned data structure; more complex queries will be shown later in the 
vignette.

```{r}
results <- ref_get("SearchCaseStudies", query=list(ID=941))
print(results)
```

You will note that the function returns a nested tibble - that is a tibble with
other data frames inside it. This means that you can interrogate the tibble as
per usual:

```{r}
cat(results[[1, "CaseStudyId"]])
cat(results[[1, "Title"]])
cat(strtrim(results[[1, "ImpactSummary"]], width = 200), "<truncated>")
cat(strtrim(results[[1, "ImpactDetails"]], width = 200), "<truncated>")
cat(results[[1, "Institution"]])
```

You can also interrogate the nested fields the same way, and even subset them:

```{r}
print(results[[1, "Country"]])
print(results[[1, "Institutions"]])
print(results[[1, "Institutions"]][,c("UKPRN", "InstitutionName")])
```

> _In the opinion of the package author, the nested tibble offers many
advantages over other data representations - it is a relatively straight-forward
exercise to transform the data into a set of wide or narrow tables if required._

Returning a single case study based on the ID is obviously a niche use-case, so
there are some other ways to search the database. But before getting to those, 
it is worth pointing out that you can select multiple case studies in a single
query:

```{r}
results <- ref_get("SearchCaseStudies", query=list(ID=c(941, 942, 1014)))
print(results)
```

The ID parameter above is an exclusive parameter - if you provide one or more 
IDs then the function will print a warning to the console, and remove all 
parameters except for the IDs. This is based on the API's documented 
limitations.

The other parameters can all be combined for searching. Those parameters are:

* **UKPRN** - This is a code referencing an institution, and comes from the 
  ListInstitutions method below. Takes a single UKPRN.
* **UoA** - This is a code referencing a Unit of Assessment, and comes from the
  ListUnitsOfAssessment method below. Takes a single ID.
* **tags** - This is one or more codes referencing tags from the ListTagValues 
  method. The tags are separated into 13 different TagTypes, which are detailed
  below. When multiple tags are provided to the search method, it will only 
  return rows which contain both tags.
* **phrase** - You can search the database using a text query. The query must
  conform to Lucene search query syntax.
  
Some examples are shown below.

```{r}
results <- ref_get("SearchCaseStudies", query=list(UKPRN = 10007777))
dim(results)
results <- ref_get("SearchCaseStudies", query=list(UoA = 5))
dim(results)
results <- ref_get("SearchCaseStudies", query=list(tags = c(11280, 5085)))
dim(results)
results <- ref_get("SearchCaseStudies", query=list(phrase = "hello"))
dim(results)
results <- ref_get("SearchCaseStudies", query=list(UKPRN = 10007146,
                                                   UoA   = 3))
dim(results)
```

Unfortunately, the API method requires at least one search parameter, which 
makes it more difficult to download the entire dataset. A short script for this
purpose is included at the end of this vignette.

Useful values for the UKPRN, UoA and tags parameters can be found by querying 
the other 4 API methods - the phrase parameter is the only parameter which can 
be used in isolation. Each of the 4 other API methods are outlined below.

## ListInstitutions

This method lists all of the institutions which are included in the REF Impact
Case Studies database, and the UKPRN column in the resuling tibble can be used
as a query parameter

```{r}
institutions <- ref_get("ListInstitutions")
print(institutions)
```

## ListTagTypes and ListTagValues

These methods provide tags which can be used as search parameters in the 
SearchCaseStudies method. The ListTagTypes method returns the types of tags 
available:

```{r}
tag_types <- ref_get("ListTagTypes")
print(tag_types)
```

These tag types can then be used as an argument to the ListTagValues method, to
get all tags for each type:

```{r}
tag_values_5 <- ref_get("ListTagValues", tag_type = 5)
print(tag_values_5)
```

This can take some time to iterate through, so the full table is bundled with 
this package. You can access it via `ref_tags`:

```{r}
print(ref_tags)
```

## ListUnitsOfAssessment

This method lists all of the units of assessment which the Impact
Case Studies can be assessed against. The tibble also includes an ID column
which can be used when querying the SearchCaseStudies method.

```{r}
UoAs <- ref_get("ListUnitsOfAssessment")
print(UoAs)
```

## Extracting the entire dataset

As alluded to above, the API cannot be searched without parameters, which means
that downloading the entire dataset is not a simple task. The code below
can be used to extract all records from the database.

```{r, eval = F}
uoa_table <- ref_get("ListUnitsOfAssessment")
uoa_list <- uoa_table$ID

ref_corpus <- vector(length = length(uoa_list), mode = "list")

for (i in seq_along(uoa_list)) {
  message("Retrieving data for UoA ", uoa_list[i])
  ref_corpus[[i]] <- ref_get("SearchCaseStudies", query = list(UoA = uoa_list[i]))
}

output <- do.call(rbind, ref_corpus)
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/x_search_case_studies.R
\name{get_case_studies}
\alias{get_case_studies}
\title{Search and download case studies}
\usage{
get_case_studies(ID = NULL, UKPRN = NULL, UoA = NULL, tags = NULL,
  phrase = NULL)
}
\arguments{
\item{ID}{integer, can return multiple case studies if provided as a numeric
vector. If this argument is provided then all other arguments will be
ignored.}

\item{UKPRN}{integer, filter for a specific institution using its UKPRN.}

\item{UoA}{integer, filter for a specific Unit of Assessment.}

\item{tags}{integer, filter for specific tag IDs. Can combine multiple tags
if provided as a numeric vector. Multiple tags are combined with logical
AND.}

\item{phrase}{string, filter for a specific phrase. Currently only supports
single words.}
}
\value{
Returns a data_frame (from the \code{tibble} package).
}
\description{
DEPRECATED - USE ref_get. This function uses the \code{SearchCaseStudies}
method from the database API. The method requires at least one filtering
parameter, which means you need to provide at least one argument to this
function.
}
\details{
This function returns a data_frame (from the \code{tibble} package) as it
deals nicely with the nested data structures provided by the API. See the
example code below for a demonstration of this behaviour.
}
\examples{
\dontrun{
studies <- get_case_studies(ID = c(27,29))
studies
studies[1,"Continent"] # [] gives list element
studies[[1,"Continent"]] # [[]] gives contents
}


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/refimpact.R
\docType{package}
\name{refimpact}
\alias{refimpact}
\alias{refimpact-package}
\title{refimpact: API wrapper for the UK REF 2014 Impact Case Studies Database}
\description{
This package provides wrapper functions around the UK Research Excellence
Framework 2014 Impact Case Studies Database API
(\url{http://impact.ref.ac.uk/}). The database contains relevant publication
and research metadata about each case study as well as several paragraphs of
text from the case study submissions. Case studies in the database are
licenced under a CC-BY 4.0 licence
(\url{http://creativecommons.org/licenses/by/4.0/legalcode}).
}
\details{
To get started, see \code{\link{ref_get}} and read the vignette using
\code{vignette('refimpact')}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/x_get_tag_types.R
\name{get_tag_types}
\alias{get_tag_types}
\title{List Tag Types}
\usage{
get_tag_types()
}
\value{
Returns a data_frame (from the \code{tibble} package).
}
\description{
DEPRECATED - USE ref_get. This function uses the \code{ListTagTypes} method
from the database API and returns a list of Tag Types and associated
metadata.
}
\examples{
\dontrun{
get_tag_types()
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/refimpact.R
\docType{data}
\name{ref_tags}
\alias{ref_tags}
\title{Table of all REF Impact Case Study tags.}
\format{A data frame with 9,400 rows and 4 variables:
\describe{
  \item{ID}{integer, the Tag ID for use in the SearchCaseStudies method}
  \item{Name}{string, the name of the tag}
  \item{TypeID}{integer, the Tag Type ID used in the ListTagValues method}
  \item{TagType}{string, the Tag Type for each tag}
}}
\source{
\url{http://impact.ref.ac.uk/}
}
\usage{
ref_tags
}
\description{
This table contains the complete set of REF Impact Case Study tags, saving
the user the effort of iteratively querying the API for each tag type.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/x_get_institutions.R
\name{get_institutions}
\alias{get_institutions}
\title{List Institutions}
\usage{
get_institutions()
}
\value{
Returns a data_frame (from the \code{tibble} package).
}
\description{
DEPRECATED - USE ref_get. This function uses the \code{ListInstitutions}
method from the database API and returns a list of institutions and
associated metadata, including the UKPRN.
}
\examples{
\dontrun{
get_institutions()
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/x_get_units_of_assessment.R
\name{get_units_of_assessment}
\alias{get_units_of_assessment}
\title{List Units of Assessment}
\usage{
get_units_of_assessment()
}
\value{
Returns a data_frame (from the \code{tibble} package).
}
\description{
DEPRECATED - USE ref_get. This function uses the \code{ListUnitsOfAssessment}
method from the database API and returns a list of Units of Assessment and
associated metadata.
}
\examples{
\dontrun{
get_units_of_assessment()
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/x_get_tag_values.R
\name{get_tag_values}
\alias{get_tag_values}
\title{List Tag Values}
\usage{
get_tag_values(tag_type)
}
\arguments{
\item{tag_type}{numeric, a valid tag type ID. Use \code{get_tag_types()} to find valid
tag types.}
}
\value{
Returns a data_frame (from the \code{tibble} package).
}
\description{
DEPRECATED - USE ref_get. This function uses the \code{ListTagValues} method
from the database API and returns a list of Tag Values for the supplied tag
type.
}
\examples{
\dontrun{
get_tag_values(3) # Tag Type 3 (Subject)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ref_get.R
\name{ref_get}
\alias{ref_get}
\title{Call the REF Impact Case Studies API}
\usage{
ref_get(api_method, tag_type = NULL, query = NULL)
}
\arguments{
\item{api_method}{text, the API method you wish to call. Valid methods
are summarised below, and documented on the REF Impact Case Studies website
linked above, as well as in the vignette.}

\item{tag_type}{integer, for ListTagValues method only. This is the ID of the
tag type you wish to retrieve. See example usage below.}

\item{query}{list, search parameters for use with the SearchCaseStudies
method. See example usage below.}
}
\value{
Returns a \code{\link[tibble]{tibble}} with nested data frames. To
  access the nested data frames, subset the tibble using the [[]] syntax. For
  more information, see the vignette.
}
\description{
This function calls the REF Impact Case Studies API, and returns the dataset
as a tibble. See the vignette for more details about how to use this
function.
}
\details{
Details about the API can be found at
\url{http://impact.ref.ac.uk/CaseStudies/APIhelp.aspx}.
}
\section{Valid API methods}{

\itemize{
  \item ListInstitutions (no arguments)
  \item ListTagTypes (no arguments)
  \item ListTagValues (tag_type is a compulsory argument)
  \item ListUnitsOfAssessment (no arguments)
  \item SearchCaseStudies (query is a compulsory argument - see below)
}
}

\section{SearchCaseStudies query argument}{

This argument is used to pass search parameters through to the API. These
parameters are passed as a named list, and you must provide at least one
parameter for this method. There are 5 parameters:
\itemize{
\item ID - Takes a single ID or a vector of IDs. If you use this parameter
you cannot use any of the other 4 parameters.
\item UKPRN (UK Provider Reference Number) - takes a single UKPRN. You can
get a list of valid values using the ListInstitutions method.
\item UoA - This is a code referencing a Unit of Assessment, and you can get
a list of valid values from the ListUnitsOfAssessment method. Takes a single
UoA.
\item tags - This is one or more codes referencing tags from the
ListTagValues method. When multiple tags are provided to the search method,
it will only return rows which contain both tags. To help you discover tags
that you can use here, you can look at the ref_tags dataset (bundled with
this package)
\item phrase - You can search the database using a text query. The query must
conform to Lucene search query syntax.
}
For more information about how to use these parameters, see the vignette.
}

\examples{
\donttest{
institutions <- ref_get("ListInstitutions")
units_of_assessment <- ref_get("ListUnitsOfAssessment")
tag_types <- ref_get("ListTagTypes")
tag_type_5 <- ref_get("ListTagValues", 5L)
ref_get("SearchCaseStudies", query = list(ID     = c(27121,1698)))
ref_get("SearchCaseStudies", query = list(UKPRN  = 10007777))
ref_get("SearchCaseStudies", query = list(UoA    = 5))
ref_get("SearchCaseStudies", query = list(tags   = c(11280, 5085)))
ref_get("SearchCaseStudies", query = list(phrase = "hello"))
ref_get("SearchCaseStudies", query = list(UKPRN  = 10007146, UoA = 3))
}

}
