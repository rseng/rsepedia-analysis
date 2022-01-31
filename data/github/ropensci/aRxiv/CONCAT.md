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
# aRxiv

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](
https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/aRxiv)](https://cran.r-project.org/package=aRxiv)

## R interface to arXiv

[arXiv](https://arxiv.org) is a repository of electronic preprints for
computer science, mathematics, physics, quantitative biology,
quantitative finance, and statistics. The
[aRxiv](https://github.com/ropensci/aRxiv) package is an R interface to
the [arXiv API](https://arxiv.org/help/api/index).

Note that the arXiv API _does not_ require an API key.


## Package Status and Installation

[![R-CMD-check](https://github.com/ropensci/aRxiv/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/aRxiv/actions)
[![codecov](https://codecov.io/gh/ropensci/aRxiv/branch/master/graph/badge.svg)](https://app.codecov.io/gh/ropensci/aRxiv)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/aRxiv?color=blue)](https://github.com/r-hub/cranlogs.app)

__Installation instructions__
__Stable Version__

You can install the package via [CRAN](https://cran.r-project.org):

```r
install.packages("aRxiv")
```

__Development Version__

Or use `devtools::install_github()` to get the (more recent) version
at [GitHub](https://github.com/rOpenSci/aRxiv):

```r
install.packages("devtools")
library(devtools)
install_github("ropensci/aRxiv")
```

## Usage
### Basic usage

The main function is `arxiv_search()`. Here's an example of its use:

```r
library(aRxiv)
z <- arxiv_search(query = 'au:"Peter Hall" AND cat:stat*', limit=50)
str(z)
```


### Tutorial

An aRxiv tutorial is available at the rOpenSci website, [here](https://docs.ropensci.org/aRxiv/articles/aRxiv.html).

To view the tutorial from R, use:

```r
vignette("aRxiv", "aRxiv")
```


### Links

* [arXiv](https://arxiv.org)
* [arXiv API](https://arxiv.org/help/api/index)
* [arXiv API user manual](https://arxiv.org/help/api/user-manual)
* [Bulk data access to arXiv](https://arxiv.org/help/bulk_data)
* [Bulk data access to arXiv metadata via OAI-PMH](https://arxiv.org/help/oa/index)
* [Bulk data access to arXiv PDFs and source docs](https://arxiv.org/help/bulk_data_s3)


### License

Licensed under the [MIT license](https://cran.r-project.org/web/licenses/MIT). ([More information here](https://en.wikipedia.org/wiki/MIT_License).)

---

This package is part of a richer suite called [fulltext](https://github.com/ropensci/fulltext), along with several other packages, that provides the ability to search for and retrieve full text of open access scholarly articles. We recommend using `fulltext` as the primary R interface to `arXiv` unless your needs are limited to this single source.

---

## Citation

Get citation information for `aRxiv` in R by running: `citation(package = 'aRxiv')`

## Code of Conduct

Please note that this project is released with a [Contributor Code of Conduct](https://github.com/ropensci/aRxiv/blob/master/CONDUCT.md).
By participating in this project you agree to abide by its terms.



[![ropensci footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
## Script for [aRxiv](http://github.com/ropensci/aRxiv) package

[`grab_api_manual_tables.R`](http://github.com/ropensci/aRxiv/tree/master/inst/scripts/grab_api_manual_tables.R)
&ndash; R script to grab tables from the
[arXiv API user guide](http://arxiv.org/help/api/index).

- search terms (`query_prefixes`)
- subject classifications (`arxiv_cats`)

The script creates datasets for the package that contain the body of the tables.

To access the resulting datasets, do the following:

```r
library(aRxiv)
data(query_prefixes)
data(arxiv_cats)
```
---
title: aRxiv tutorial
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{aRxiv tutorial}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8](inputenc)
---

[arXiv](https://arxiv.org) is a repository of electronic preprints for
computer science, mathematics, physics, quantitative biology,
quantitative finance, and statistics. The
[aRxiv package](https://github.com/ropensci/aRxiv) provides an
[R](https://www.r-project.org) interface to the
[arXiv API](https://arxiv.org/help/api/index).

Note that the arXiv API _does not_ require an API key.

```{r be_graceful_on_cran, include=FALSE}
# if on cran, avoid errors when connection problems
on_cran <- Sys.getenv("NOT_CRAN")!="true"
if(on_cran) {
    aRxiv:::set_arxiv_timeout(1)
    aRxiv:::set_message_on_timeout(TRUE)
}
```

```{r change_aRxiv_delay_option, include=FALSE}
options(aRxiv_delay=0.5)
```

## Installation

You can install the [aRxiv package](https://github.com/rOpenSci/aRxiv)
via [CRAN](https://cran.r-project.org):

```{r install_from_cran, eval=FALSE}
install.packages("aRxiv")
```

Or use `devtools::install_github()` to get the (possibly more recent) version
at [GitHub](https://github.com/rOpenSci/aRxiv):

```{r install_pkgs, eval=FALSE}
install.packages("devtools")
library(devtools)
install_github("ropensci/aRxiv")
```

## Basic use

Use `arxiv_search()` to search [arXiv](https://arxiv.org),
`arxiv_count()` to get a simple count of manuscripts matching a
query, and `arxiv_open()` to open the abstract pages for a set of
results from `arxiv_search()`.

We'll get to the details in a moment. For now, let's look at a few
examples.

Suppose we wanted to identify all arXiv manuscripts with &ldquo;`Peter
Hall`&rdquo; as an author. It is best to first get a count, so that we have
a sense of how many records the search will return. (Peter Hall was
&ldquo;[among the world's most prolific and highly cited authors in both probability and statistics](https://en.wikipedia.org/wiki/Peter_Gavin_Hall).&rdquo;)
We first use `library()` to load the aRxiv package and then
`arxiv_count()` to get the count.

```{r arxiv_count}
library(aRxiv)
arxiv_count('au:"Peter Hall"')
```

The `au:` part indicates to search the author field; we use double
quotes to search for a _phrase_.

To obtain the actual records matching the query, use `arxiv_search()`.

```{r arxiv_search}
rec <- arxiv_search('au:"Peter Hall"')
nrow(rec)
```

The default is to grab no more than 10 records; this limit can be
changed with the `limit` argument. But note that the arXiv API will
not let you download more than 50,000 or so records, and even in that
case it's best to do so in batches; more on this below.

Also note that the result of `arxiv_search()` has an attribute
`"total_results"` containing the total count of search results; this
is the same as what `arxiv_count()` provides.

```{r arxiv_search_attr}
attr(rec, "total_results")
```

The following will get us all `r arxiv_count('au:"Peter Hall"')`
records.

```{r arxiv_search_limit100}
rec <- arxiv_search('au:"Peter Hall"', limit=100)
nrow(rec)
```

`arxiv_search()` returns a data frame with each row being a single
manuscript. The columns are the different fields (e.g., `authors`, `title`,
`abstract`, etc.). Fields like `authors` that
contain multiple items will be a single character string with the
multiple items separated by a vertical bar (`|`).

We might be interested in a more restrictive search, such as for Peter
Hall's arXiv manuscripts that have `deconvolution` in the title. We
use `ti:` to search the title field, and combine the two with `AND`.

```{r arxiv_search_deconvolution}
deconv <- arxiv_search('au:"Peter Hall" AND ti:deconvolution')
nrow(deconv)
```

Let's display just the authors and title for the results.

```{r authors_title}
deconv[, c('title', 'authors')]
```

We can open the abstract pages for these `r nrow(deconv)` manuscripts
using `arxiv_open()`. It takes, as input, the output of
`arxiv_search()`.

```{r arxiv_open, eval=FALSE}
arxiv_open(deconv)
```


## Forming queries

The two basic arguments to `arxiv_count()` and `arxiv_search()` are
`query`, a
character string representing the search, and `id_list`, a list of
[arXiv manuscript identifiers](https://arxiv.org/help/arxiv_identifier).

- If only `query` is provided, manuscripts matching that query are
  returned.
- If only `id_list` is provided, manuscripts in the list are
  returned.
- If both are provided, manuscripts in `id_list` that match `query`
  will be returned.

`query` may be a single character string or a vector of character
strings. If it is a vector, the elements are pasted together with
`AND`.

`id_list` may be a vector of character strings or a single
comma-separated character string.

### Search terms

Generally, one would ignore `id_list` and focus on forming the `query`
argument. The aRxiv package includes a dataset `query_terms` that
lists the terms (like `au`) that you can use.

```{r query_terms}
query_terms
```

Use a colon (`:`) to separate the query term from the actual query.
Multiple queries can be combined with `AND`, `OR`, and `ANDNOT`. The
default is `OR`.

```{r illustrate_AND}
arxiv_count('au:Peter au:Hall')
arxiv_count('au:Peter OR au:Hall')
arxiv_count('au:Peter AND au:Hall')
arxiv_count('au:Hall ANDNOT au:Peter')
```

It appears that in the author field (and many other fields) you must
search full words, and that wild cards are not allowed.

```{r illustrate_wildcard}
arxiv_count('au:P* AND au:Hall')
arxiv_count('au:P AND au:Hall')
arxiv_count('au:"P Hall"')
```

### Subject classifications

arXiv has a set of `r nrow(arxiv_cats)` subject classifications,
searchable with the prefix `cat:`. The aRxiv package contains a
dataset `arxiv_cats` containing the abbreviations and descriptions.
Here are the statistics categories.

```{r arxiv_cats}
arxiv_cats[grep('^stat', arxiv_cats$abbreviation),]
```

To search these categories, you need to include either the full term
or use the `*` wildcard.

```{r search_cats}
arxiv_count('cat:stat')
arxiv_count('cat:stat.AP')
arxiv_count('cat:stat*')
```

### Dates and ranges of dates

The terms `submittedDate` (date/time of first submission) and
`lastUpdatedDate` (date/time of last revision) are particularly
useful for limiting a search with _many_ results, so that you may
combine multiple searches together, each within some window of time,
to get the full results.

The date/time information is of the form `YYYYMMDDHHMMSS`, for
example `20071018122534` for `2007-10-18 12:25:34`. You can use `*`
for a wildcard for the times. For example, to get all manuscripts
with initial submission on 2007-10-18:

```{r wildcard_times}
arxiv_count('submittedDate:20071018*')
```

But you can't use the wildcard within the _dates_.

```{r wildcard_date}
arxiv_count('submittedDate:2007*')
```

To get a count of all manuscripts with original submission in 2007,
use a date range, like `[from_date TO to_date]`. (If you give a partial
date, it's treated as the earliest date/time that matches, and the
range appears to be up to but not including the second date/time.)

```{r daterange}
arxiv_count('submittedDate:[2007 TO 2008]')
```

## Search results

The output of `arxiv_search()` is a data frame with the following
columns.

```{r arxiv_search_result}
res <- arxiv_search('au:"Peter Hall"')
names(res)
```

The columns are described in the help file for `arxiv_search()`. Try
`?arxiv_search`.

A few short notes:

- Each field is a single character string. `authors`, `link_doi`, and
  `categories` may contain multiple items, separated by a vertical bar
  (`|`).
- Missing entries will have an empty character string (`""`).
- The `categories` column may contain not just the aRxiv categories
  (e.g., `stat.AP`) but also codes for the
  [Mathematical Subject Classification (MSC)](https://mathscinet.ams.org/mathscinet/msc/msc2010.html) (e.g., 14J60)
  and the
  [ACM Computing Classification System](https://www.acm.org/about/class/1998/)
  (e.g., F.2.2). These are not searchable with `cat:` but are
  searchable with a general search.

```{r search_msc}
arxiv_count("cat:14J60")
arxiv_count("14J60")
```


## Sorting results

The `arxiv_search()` function has two arguments for sorting the results,
`sort_by` (taking values `"submitted"`, `"updated"`, or
`"relevance"`) and `ascending` (`TRUE` or `FALSE`). If `id_list` is
provided, these sorting arguments are ignored and the results are
presented according to the order in `id_list`.

Here's an example, to sort the results by the date the manuscripts
were last updated, in descending order.

```{r sortby_example}
res <- arxiv_search('au:"Peter Hall" AND ti:deconvolution',
                    sort_by="updated", ascending=FALSE)
res$updated
```


## Technical details

### Metadata limitations

The [arXiv metadata](https://arxiv.org/help/prep) has a number of
limitations, the key issue being that it is author-supplied and so not
necessarily consistent between records.

Authors' names may vary between records (e.g., Peter Hall vs. Peter G.
Hall vs. Peter Gavin Hall vs. P Hall). Further, arXiv provides no ability to
distinguish multiple individuals with the same name (c.f.,
[ORCID](https://orcid.org)).

Authors' institutional affiliations are mostly missing.  The arXiv
submission form does not include an affiliation field; affiliations
are entered within the author field, in parentheses.  The
[metadata instructions](https://arxiv.org/help/prep#author) may not be
widely read.

There are no key words; you are stuck with searching the free text in
the titles and abstracts.

Subject classifications are provided by the authors and may be
incomplete or inappropriate.


### Limit time between search requests

Care should be taken to avoid multiple requests to the arXiv API in a
short period of time. The
[arXiv API user manual](https://arxiv.org/help/api/user-manual) states:

> In cases where the API needs to be called multiple times in a row,
> we encourage you to play nice and incorporate a 3 second delay in
> your code.

The aRxiv package institutes a delay between requests, with the time
period for the delay configurable with the R option
`"aRxiv_delay"` (in seconds). The default is 3 seconds.

To reduce the delay to 1 second, use:

```{r aRxiv_delay, eval=FALSE}
options(aRxiv_delay=1)
```

**Don't** do searches in parallel (e.g., via the parallel
package). You may be locked out from the arXiv API.


### Limit number of items returned

The arXiv API returns only complete records (including the entire
abstracts); searches returning large numbers of records can be very
slow.

It's best to use `arxiv_count()` before `arxiv_search()`, so that you have
a sense of how many records you will receive. If the count is large,
you may wish to refine your query.

arXiv has a hard limit of around 50,000 records; for a query that
matches more than 50,000 manuscripts, there is no way to receive the
full results. The simplest solution to this problem is to break the
query into smaller pieces, for example using slices of time, with a
range of dates for `submittedDate` or `lastUpdatedDate`.

The `limit` argument to `arxiv_search()` (with default `limit=10`)
limits the number of records to be returned. If you wish to receive
more than 10 records, you must specify a larger limit (e.g., `limit=100`).

To avoid accidental searches that may return a very large number of
records, `arxiv_search()` uses an R option, `aRxiv_toomany` (with a
default of 15,000), and refuses to attempt a search that will return
results above that limit.

### Make requests in batches

Even for searches that return a moderate number of records (say
2,000), it may be best to make the requests in batches: Use a smaller
value for the `limit` argument (say 100), and make multiple requests
with different offsets, indicated with the `start` argument, for the
initial record to return.

This is done automatically with the `batchsize` argument to
`arxiv_search()`.  A search is split into multiple calls, with no more
than `batchsize` records to be returned by each, and then the results
are combined.


## License and bugs

- License:
  [MIT](https://github.com/ropensci/aRxiv/blob/master/LICENSE)
- Report bugs or suggestions improvements by [submitting an issue](https://github.com/ropensci/aRxiv/issues) to
  [our GitHub repository for aRxiv](https://github.com/ropensci/aRxiv).


```{r reset_to_defaults, include=FALSE}
if(on_cran) {
    aRxiv:::set_arxiv_timeout(30)
    aRxiv:::set_message_on_timeout(FALSE)
}
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/query_terms-data.R
\docType{data}
\name{query_terms}
\alias{query_terms}
\title{arXiv query field terms}
\format{
A data frame with two columns: the \code{term} and corresponding
\code{description}.
}
\source{
\url{https://arxiv.org/help/api/user-manual}
}
\usage{
data(query_terms)
}
\description{
Possible terms that correspond to different fields in arXiv searches.
}
\examples{
query_terms
}
\author{
Karl W Broman
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/arxiv_search.R
\name{arxiv_search}
\alias{arxiv_search}
\title{The main search function for aRxiv}
\usage{
arxiv_search(
  query = NULL,
  id_list = NULL,
  start = 0,
  limit = 10,
  sort_by = c("submitted", "updated", "relevance"),
  ascending = TRUE,
  batchsize = 100,
  force = FALSE,
  output_format = c("data.frame", "list"),
  sep = "|"
)
}
\arguments{
\item{query}{Search pattern as a string; a vector of such strings
also allowed, in which case the elements are combined with \code{AND}.}

\item{id_list}{arXiv doc IDs, as comma-delimited string or a vector
of such strings}

\item{start}{An offset for the start of search}

\item{limit}{Maximum number of records to return.}

\item{sort_by}{How to sort the results (ignored if \code{id_list} is
provided)}

\item{ascending}{If TRUE, sort in ascending order; else descending
(ignored if \code{id_list} is provided)}

\item{batchsize}{Maximum number of records to request at one time}

\item{force}{If TRUE, force search request even if it seems extreme}

\item{output_format}{Indicates whether output should be a data frame or a list.}

\item{sep}{String to use to separate multiple authors,
affiliations, DOI links, and categories, in the case that
\code{output_format="data.frame"}.}
}
\value{
If \code{output_format="data.frame"}, the result is a data
frame with each row being a manuscript and columns being the
various fields.

If \code{output_format="list"}, the result is a list parsed from
the XML output of the search, closer to the raw output from arXiv.

The data frame format has the following columns.
\tabular{rll}{
[,1] \tab id               \tab arXiv ID \cr
[,2] \tab submitted        \tab date first submitted \cr
[,3] \tab updated          \tab date last updated \cr
[,4] \tab title            \tab manuscript title \cr
[,5] \tab summary          \tab abstract \cr
[,6] \tab authors          \tab author names \cr
[,7] \tab affiliations     \tab author affiliations \cr
[,8] \tab link_abstract    \tab hyperlink to abstract \cr
[,9] \tab link_pdf         \tab hyperlink to pdf \cr
[,10] \tab link_doi         \tab hyperlink to DOI \cr
[,11] \tab comment          \tab authors' comment \cr
[,12] \tab journal_ref      \tab journal reference \cr
[,13] \tab doi              \tab published DOI \cr
[,14] \tab primary_category \tab primary category \cr
[,15] \tab categories       \tab all categories \cr
}

The contents are all strings; missing values are empty strings (\code{""}).

The columns \code{authors}, \code{affiliations}, \code{link_doi},
and \code{categories} may have multiple entries separated by
\code{sep} (by default, \code{"|"}).

The result includes an attribute \code{"search_info"} that includes
information about the details of the search parameters, including
the time at which it was completed. Another attribute
\code{"total_results"} is the total number of records that match
the query.
}
\description{
Allows for progammatic searching of the arXiv pre-print repository.
}
\examples{
\dontshow{old_delay <- getOption("aRxiv_delay")
          options(aRxiv_delay=1)}
\donttest{
# search for author Peter Hall with deconvolution in title
z <- arxiv_search(query = 'au:"Peter Hall" AND ti:deconvolution', limit=2)
attr(z, "total_results") # total no. records matching query
z$title

# search for a set of documents by arxiv identifiers
z <- arxiv_search(id_list = c("0710.3491v1", "0804.0713v1", "1003.0315v1"))
# can also use a comma-separated string
z <- arxiv_search(id_list = "0710.3491v1,0804.0713v1,1003.0315v1")
# Journal references, if available
z$journal_ref

# search for a range of dates (in this case, one day)
z <- arxiv_search("submittedDate:[199701010000 TO 199701012400]", limit=2)
}
\dontshow{options(aRxiv_delay=old_delay)}

}
\seealso{
\code{\link[=arxiv_count]{arxiv_count()}}, \code{\link[=arxiv_open]{arxiv_open()}},
\code{\link[=query_terms]{query_terms()}}, \code{\link[=arxiv_cats]{arxiv_cats()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/can_arxiv_connect.R
\name{can_arxiv_connect}
\alias{can_arxiv_connect}
\title{Check for connection to arXiv API}
\usage{
can_arxiv_connect(max_time = 5)
}
\arguments{
\item{max_time}{Maximum wait time in seconds}
}
\value{
Returns TRUE if connection is established and FALSE
otherwise.
}
\description{
Check for connection to arXiv API
}
\examples{
\donttest{
can_arxiv_connect(2)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/arxiv_count.R
\name{arxiv_count}
\alias{arxiv_count}
\title{Count number of results for a given search}
\usage{
arxiv_count(query = NULL, id_list = NULL)
}
\arguments{
\item{query}{Search pattern as a string; a vector of such strings is
also allowed, in which case the elements are combined with \code{AND}.}

\item{id_list}{arXiv doc IDs, as comma-delimited string or a vector
of such strings}
}
\value{
Number of results (integer). An attribute
\code{"search_info"} contains information about the search
parameters and the time at which it was performed.
}
\description{
Count the number of results for a given search. Useful to check
before attempting to pull down a very large number of records.
}
\examples{
\dontshow{old_delay <- getOption("aRxiv_delay")
          options(aRxiv_delay=1)}
\donttest{
# count papers in category stat.AP (applied statistics)
arxiv_count(query = "cat:stat.AP")

# count papers by Peter Hall in any stat category
arxiv_count(query = 'au:"Peter Hall" AND cat:stat*')

# count papers for a range of dates
#    here, everything in 2013
arxiv_count("submittedDate:[2013 TO 2014]")
}
\dontshow{options(aRxiv_delay=old_delay)}

}
\seealso{
\code{\link[=arxiv_search]{arxiv_search()}}, \code{\link[=query_terms]{query_terms()}},
\code{\link[=arxiv_cats]{arxiv_cats()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/arxiv_cats-data.R
\docType{data}
\name{arxiv_cats}
\alias{arxiv_cats}
\title{arXiv subject classifications}
\format{
A data frame with two columns: the abbreviations of the
subject classifications (\code{abbreviation}) and the corresponding
description (\code{description}).
}
\source{
\url{https://arxiv.org/help/api/user-manual}
}
\usage{
data(arxiv_cats)
}
\description{
arXiv subject classifications: their abbreviations and
corresponding descriptions.
}
\examples{
arxiv_cats
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/arxiv_open.R
\name{arxiv_open}
\alias{arxiv_open}
\title{Open abstract for results of arXiv search}
\usage{
arxiv_open(search_results, limit = 20)
}
\arguments{
\item{search_results}{Data frame of search results, as returned from \code{\link[=arxiv_search]{arxiv_search()}}.}

\item{limit}{Maximum number of abstracts to open in one call.}
}
\value{
(Invisibly) Vector of character strings with URLs of
abstracts opened.
}
\description{
Open, in web browser, the abstract pages for each of set of arXiv search results.
}
\details{
There is a delay between calls to
\code{\link[utils:browseURL]{utils::browseURL()}}, with the amount taken from the R
option \code{"aRxiv_delay"} (in seconds); if missing, the default
is 3 sec.
}
\examples{
\donttest{z <- arxiv_search('au:"Peter Hall" AND ti:deconvolution')
arxiv_open(z)}

}
\seealso{
\code{\link[=arxiv_search]{arxiv_search()}}
}
