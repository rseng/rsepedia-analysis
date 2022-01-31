
<!-- README.md is generated from README.Rmd. Please edit that file -->

# internetarchive

[![Build
Status](https://travis-ci.org/ropensci/internetarchive.svg?branch=master)](https://travis-ci.org/ropensci/internetarchive)
[![codecov.io](https://codecov.io/github/ropensci/internetarchive/coverage.svg?branch=master)](https://codecov.io/github/ropensci/internetarchive?branch=master)

## Overview

This API client for the [Internet Archive](https://archive.org/) is
intended primarily for searching for items, retrieving metadata for
items, and downloading the files associated with items. The functions
can be used with the pipe operator (`%>%`) from
[magrittr](https://github.com/smbache/magrittr) and the data
manipulation verbs in [dplyr](https://github.com/hadley/dplyr) to create
pipelines from searching to downloading. For the full details of what is
possible with the Internet Archive API, see their [advanced search
help](https://archive.org/advancedsearch.php).

## Installation

Install this package from CRAN:

``` r
install.packages("internetarchive")
```

Or, install the [development
version](https://github.com/ropensci/internetarchive) from GitHub with
[devtools](http://cran.rstudio.org/web/packages/devtools/).

``` r
# install.packages("devtools")
devtools::install_github("ropensci/internetarchive", build_vignettes = TRUE)
```

Then load the package. We will also use
[dplyr](https://github.com/hadley/dplyr) for manipulating the retrieved
data.

``` r
library("internetarchive")
library("dplyr", warn.conflicts = FALSE)
```

## Basic search and browse

The simplest way to search the Internet Archive is to use a keyword
search. The following function searches for these keywords in the most
important metadata fields, and returns a list of item identifiers.

``` r
ia_keyword_search("isaac hecker")
#> 31 total items found. This query requested 5 results.
#> [1] "americanexperien00fari"    "fatherhecker00sedggoog"   
#> [3] "fatherhecker01sedg"        "abitunpublished00heckgoog"
#> [5] "TheLifeOfFatherHecker"
```

You can pass an item identifier to the `ia_browse()` function to open an
item in your browser. If you pass this function multiple identifiers, it
will open only the first one.

``` r
ia_browse("TheLifeOfFatherHecker")
```

## Advanced search

Usually it is more useful to perform an advanced search. You can
construct an advanced search as a named character vector, where the
names correspond to the fields. The following search, for instance,
looks for items published by the American Tract Society in 1864. Run the
function `ia_list_fields()` to see the list of accepted metadata fields.

``` r
ats_query <- c("publisher" = "american tract society", "year" = "1864")
ia_search(ats_query, num_results = 20)
#> 13 total items found. This query requested 20 results.
#>  [1] "huguenotsfrance00martgoog" "missionsmartyrsi00bost"   
#>  [3] "vitalgodlinessa00plumgoog" "liliantaleofthre00lili"   
#>  [5] "littlewillietrue00amer"    "ourvillageinwart00mart"   
#>  [7] "vitalgodlinesstrws00plum"  "vitalgodlinesstr00plum"   
#>  [9] "songsofzionenlar00amer"    "ilvertonrectoryo00mart"   
#> [11] "colorbearerfranc01amer"    "sketcheseloquen00wategoog"
#> [13] "sketchesofeloque00wate"
```

You can change the number of items returned by the search using the
`num_results =` argument, and you can request subsequent pages of
results with the `page =` argument.

Notice that `ia_search()` and `ia_keyword_search()` both return a
character vector of identifiers, so both can be used in the same way at
the beginning of a pipeline.

### Dates

To search by a date range, use the `date` field and the years (or
[ISO 8601 dates](http://en.wikipedia.org/wiki/ISO_8601)) separated by
`TO`. Here we search for publications by the American Tract Society in
the
1840s.

``` r
ia_search(c("publisher" = "american tract society", date = "1840 TO 1850"))
#> 104 total items found. This query requested 5 results.
#> [1] "historyreformat09aubgoog"  "scripturebiogra00hookgoog"
#> [3] "historyreformat22aubgoog"  "memoirmrssarahl00hookgoog"
#> [5] "circulationandc00socigoog"
```

## Getting item metadata and files

Once you have retrieved a list of items, you can retrieve their metadata
and the list of files associated with the items.

To get a single item’s metadata, you can pass its identifier to the
`ia_get_items()` function.

``` r
hecker <- ia_get_items("TheLifeOfFatherHecker")
```

The result is a list where the names of items in the list are the item
identifiers, and the rest of the list is the metadata. This nested list
can be difficult to work with, so the `ia_metadata()` returns a data
frame of the metadata, and `ia_files()` returns a data frame of the
files associated with the item.

``` r
ia_metadata(hecker)
ia_files(hecker)
```

These functions can also retrieve the information for multiple items
when used in a pipeline. Here we search for all the items about Hecker,
retrieve their metadata, and turn it into a data frame. We then filter
the data frame to get only the titles.

``` r
ia_keyword_search("isaac hecker", num_results = 20) %>% 
  ia_get_items() %>% 
  ia_metadata() %>% 
  filter(field == "title") %>% 
  select(value)
```

## Downloading files

The `ia_download()` function will download all the files in a data frame
returned from `ia_files()`. This function should be used with caution,
and you should first filter the data frame to download only the files
that you wish. In the following example, we retrieve a list of all the
files associated with items published by the American Tract Society in
1864. Then we filter the list so we get only text files, then we pick
only the first text file associated with each item. Finally we download
the files to a directory we specify (in this case, a temporary
directory).

``` r
dir <- tempdir()
ia_search(ats_query) %>% 
  ia_get_items() %>% 
  ia_files() %>% 
  filter(type == "txt") %>% 
  group_by(id) %>% 
  slice(1) %>% 
  ia_download(dir = dir, overwrite = FALSE) %>% 
  glimpse()
```

Notice that `ia_download()` returns a modified version of the data frame
that was passed to it, adding a column `local_file` with the path to the
download files.

If the `overwrite =` argument is `FALSE`, then you can pass the same
data frame of files to `ia_download()` and it will download only the
files that it has not already downloaded.

-----

[![rOpenSCi
logo](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
# internetarchive 0.1.6

- Fix failing test related to change at Internet Archive

# internetarchive 0.1.5

- Update package dependencies

# internetarchive 0.1.4

- Add `rmarkdown` as a suggested package since `knitr` has changed for vignettes

# internetarchive 0.1.2

- First CRAN release
- Tests and minor bugfixes

# internetarchive 0.1.1

- First public release on GitHub with rOpenSci
- Functions for searching, downloading, metadata
This  update to the `internetarchive` package fixes a failing test as requested by CRAN maintainers. The test was failing due to a change at the Internet Archive.

## Test environments

* local OS X install: R-release
* Ubuntu 14.04, via Travis CI: R-release, R-oldrel, R-devel
* win-builder: R-devel, R-release

## R CMD check results

There were no ERRORs or WARNINGs. The only NOTE pertains to a new release of the package.
---
output: github_document
pagetitle: "An API Client for the Internet Archive"
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```

# internetarchive

[![Build Status](https://travis-ci.org/ropensci/internetarchive.svg?branch=master)](https://travis-ci.org/ropensci/internetarchive)
[![codecov.io](https://codecov.io/github/ropensci/internetarchive/coverage.svg?branch=master)](https://codecov.io/github/ropensci/internetarchive?branch=master)

## Overview

This API client for the [Internet Archive](https://archive.org/) is intended primarily for searching for items, retrieving metadata for items, and downloading the files associated with items. The functions can be used with the pipe operator (`%>%`) from [magrittr](https://github.com/smbache/magrittr) and the data manipulation verbs in [dplyr](https://github.com/hadley/dplyr) to create pipelines from searching to downloading. For the full details of what is possible with the Internet Archive API, see their [advanced search help](https://archive.org/advancedsearch.php).

## Installation

Install this package from CRAN:

```{r eval=FALSE}
install.packages("internetarchive")
```

Or, install the [development version](https://github.com/ropensci/internetarchive) from GitHub with [devtools](http://cran.rstudio.org/web/packages/devtools/).

```{r eval=FALSE}
# install.packages("devtools")
devtools::install_github("ropensci/internetarchive", build_vignettes = TRUE)
```

Then load the package. We will also use [dplyr](https://github.com/hadley/dplyr) for manipulating the retrieved data.


```{r}
library("internetarchive")
library("dplyr", warn.conflicts = FALSE)
```

## Basic search and browse

The simplest way to search the Internet Archive is to use a keyword search. The following function searches for these keywords in the most important metadata fields, and returns a list of item identifiers.

```{r}
ia_keyword_search("isaac hecker")
```

You can pass an item identifier to the `ia_browse()` function to open an item in your browser. If you pass this function multiple identifiers, it will open only the first one.

```{r, eval=FALSE}
ia_browse("TheLifeOfFatherHecker")
```

## Advanced search

Usually it is more useful to perform an advanced search. You can construct an advanced search as a named character vector, where the names correspond to the fields. The following search, for instance, looks for items published by the American Tract Society in 1864. Run the function `ia_list_fields()` to see the list of accepted metadata fields. 


```{r}
ats_query <- c("publisher" = "american tract society", "year" = "1864")
ia_search(ats_query, num_results = 20)
```

You can change the number of items returned by the search using the `num_results =` argument, and you can request subsequent pages of results with the `page =` argument.

Notice that `ia_search()` and `ia_keyword_search()` both return a character vector of identifiers, so both can be used in the same way at the beginning of a pipeline.

### Dates

To search by a date range, use the `date` field and the years (or [ISO 8601 dates](http://en.wikipedia.org/wiki/ISO_8601)) separated by `TO`. Here we search for publications by the American Tract Society in the 1840s.


```{r}
ia_search(c("publisher" = "american tract society", date = "1840 TO 1850"))
```

## Getting item metadata and files

Once you have retrieved a list of items, you can retrieve their metadata and the list of files associated with the items.

To get a single item's metadata, you can pass its identifier to the `ia_get_items()` function.


```{r, eval=FALSE}
hecker <- ia_get_items("TheLifeOfFatherHecker")
```

The result is a list where the names of items in the list are the item identifiers, and the rest of the list is the metadata. This nested list can be difficult to work with, so the `ia_metadata()` returns a data frame of the metadata, and `ia_files()` returns a data frame of the files associated with the item.


```{r, eval=FALSE}
ia_metadata(hecker)
ia_files(hecker)
```

These functions can also retrieve the information for multiple items when used in a pipeline. Here we search for all the items about Hecker, retrieve their metadata, and turn it into a data frame. We then filter the data frame to get only the titles.


```{r, eval=FALSE}
ia_keyword_search("isaac hecker", num_results = 20) %>% 
  ia_get_items() %>% 
  ia_metadata() %>% 
  filter(field == "title") %>% 
  select(value)
```

## Downloading files

The `ia_download()` function will download all the files in a data frame returned from `ia_files()`. This function should be used with caution, and you should first filter the data frame to download only the files that you wish. In the following example, we retrieve a list of all the files associated with items published by the American Tract Society in 1864. Then we filter the list so we get only text files, then we pick only the first text file associated with each item. Finally we download the files to a directory we specify (in this case, a temporary directory).


```{r eval=FALSE}
dir <- tempdir()
ia_search(ats_query) %>% 
  ia_get_items() %>% 
  ia_files() %>% 
  filter(type == "txt") %>% 
  group_by(id) %>% 
  slice(1) %>% 
  ia_download(dir = dir, overwrite = FALSE) %>% 
  glimpse()
```


Notice that `ia_download()` returns a modified version of the data frame that was passed to it, adding a column `local_file` with the path to the download files.

If the `overwrite =` argument is `FALSE`, then you can pass the same data frame of files to `ia_download()` and it will download only the files that it has not already downloaded.

---
[![rOpenSCi logo](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
---
title: "Internet Archive API Client"
author: "Lincoln Mullen"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Internet Archive API Client}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

This API client for the [Internet Archive](https://archive.org/) is intended primarily for searching for items, retrieving metadata for items, and downloading the files associated with the items. The functions can be used with the pipe operator (`%>%`) from [magrittr](https://github.com/smbache/magrittr) and the data manipulation verbs in [dplyr](https://github.com/hadley/dplyr) to create pipelines from searching to downloading. For the full details of what is possible with the Internet Archive API, see their [advanced search help](https://archive.org/advancedsearch.php).

First load the package. We will also use [dplyr](https://github.com/hadley/dplyr) for manipulating the retrieved data.

```{r message=FALSE}
library(internetarchive)
library(dplyr)
```

## Basic search and browse

The simplest way to search the Internet Archive is to use a keyword search. The following function searches for these keywords in the most important metadata fields, and returns a list of item identifiers.

```{r}
ia_keyword_search("isaac hecker")
```

You can pass an item identifier to the `ia_browse()` function to open an item in your browser. If you pass this function multiple identifiers, it will open only the first one.

```{r, eval=FALSE}
ia_browse("TheLifeOfFatherHecker")
```

## Advanced search

Usually it is more useful to perform an advanced search. You can construct an advanced search as a named character vector, where the names correspond to the fields. The following search, for instance, looks for items published by the American Tract Society in 1864. Run the function `ia_list_fields()` to see the list of accepted metadata fields. 

```{r}
ats_query <- c("publisher" = "american tract society", "year" = "1864")
ia_search(ats_query, num_results = 3)
```

You can change the number of items returned by the search using the `num_results =` argument, and you can request subsequent pages of results with the `page =` argument.

Notice that `ia_search()` and `ia_keyword_search()` both return a character vector of identifiers, so both can be used in the same way at the beginning of a pipeline.

### Dates

To search by a date range, use the `date` field and the years (or [ISO 8601 dates](http://en.wikipedia.org/wiki/ISO_8601)) separated by `TO`. Here we search for publications by the American Tract Society in the 1840s.

```{r, eval=FALSE}
ia_search(c("publisher" = "american tract society", date = "1840 TO 1850"))
```

## Getting item metadata and files

Once you have retrieved a list of items, you can retrieve their metadata and the list of files associated with the items.

To get a single item's metadata, you can pass its identifier to the `ia_get_items()` function.

```{r, eval=FALSE}
hecker <- ia_get_items("TheLifeOfFatherHecker")
```

The result is a list where the names of items in the list are the item identifiers, and the rest of the list is the metadata. This nested list can be difficult to work with, so the `ia_metadata()` returns a data frame of the metadata, and `ia_files()` returns a data frame of the files associated with the item.

```{r, eval=FALSE}
ia_metadata(hecker)
ia_files(hecker)
```

These functions can also retrieve the information for multiple items when used in a pipeline. Here we search for all the items about Hecker, retrieve their metadata, and turn it into a data frame. We then filter the data frame to get only the titles.

```{r, eval=FALSE}
ia_keyword_search("isaac hecker", num_results = 3) %>% 
  ia_get_items() %>% 
  ia_metadata() %>% 
  filter(field == "title") %>% 
  select(value)
```

## Downloading files

The `ia_download()` function will download all the files in a data frame returned from `ia_files()`. This function should be used with caution, and you should first filter the data frame to download only the files that you wish. In the following example, we retrieve a list of all the files associated with items published by the American Tract Society in 1864. Then we filter the list so we get only text files, then we pick only the first text file associated with each item. Finally we download the files to a directory we specify (in this case, a temporary directory).

```{r, eval=FALSE}
dir <- tempdir()
ia_search(ats_query, num_results = 2) %>% 
  ia_get_items() %>% 
  ia_files() %>% 
  filter(type == "txt") %>% 
  group_by(id) %>% 
  slice(1) %>% 
  ia_download(dir = dir, overwrite = FALSE) %>% 
  glimpse()
```

Notice that `ia_download()` returns a modified version of the data frame that was passed to it, adding a column `local_file` with the path to the download files.

If the `overwrite =` argument is `FALSE`, then you can pass the same data frame of files to `ia_download()` and it will download only the files that it has not already downloaded.% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internetarchive-package.R
\name{internetarchive}
\alias{internetarchive}
\title{Client for the Internet Archive API}
\description{
This client permits you to search (\link{ia_search}), retrieve item metadata
(\link{ia_metadata}) and associated files (\link{ia_files}), and download
files (\link{ia_files}) in a pipeable interface.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ia_files.R
\name{ia_files}
\alias{ia_files}
\title{Access the list of files associated with an Internet Archive item}
\usage{
ia_files(items)
}
\arguments{
\item{items}{A list describing an Internet Archive items returned from
the API.}
}
\value{
A list containing the files as a list of character vectors.
}
\description{
Access the list of files associated with an Internet Archive item
}
\examples{
\dontrun{
ats_query <- c("publisher" = "american tract society")
ids       <- ia_search(ats_query, num_results = 3)
items     <- ia_get_items(ids)
files     <- ia_files(items)
files
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ia_item_id.R
\name{ia_item_id}
\alias{ia_item_id}
\title{Access the item IDs from an Internet Archive items}
\usage{
ia_item_id(item)
}
\arguments{
\item{item}{A list describing an Internet Archive items returned from
the API. This argument is vectorized.}
}
\value{
A character vector containing the item IDs.
}
\description{
Access the item IDs from an Internet Archive items
}
\examples{
ats_query <- c("publisher" = "american tract society")
ids       <- ia_search(ats_query, num_results = 3)
items     <- ia_get_items(ids)
ia_item_id(items)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ia_metadata.R
\name{ia_metadata}
\alias{ia_metadata}
\title{Access the item metadata from an Internet Archive item}
\usage{
ia_metadata(items)
}
\arguments{
\item{items}{A list object describing an Internet Archive items returned from
the API.}
}
\value{
A data frame containing the metadata, with columns \code{id} for the
  item identifier, \code{field} for the name of the metadata field, and
  \code{value} for the metadata values.
}
\description{
Access the item metadata from an Internet Archive item
}
\examples{
ats_query <- c("publisher" = "american tract society")
ids       <- ia_search(ats_query, num_results = 3)
items     <- ia_get_items(ids)
metadata  <- ia_metadata(items)
metadata
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ia_search.R
\name{ia_search}
\alias{ia_search}
\title{Search the Internet Archive}
\usage{
ia_search(terms, num_results = 5, page = 1, print_url = FALSE,
  print_total = TRUE)
}
\arguments{
\item{terms}{A set of metadata fields and corresponding values to search.
These should take the form of a named character vector.}

\item{num_results}{The number of results to return per page.}

\item{page}{When results are paged, which page of results to return.}

\item{print_url}{Should the URL used for the query be printed as a message?}

\item{print_total}{Should the total number of results for this query be
printed as a message?}
}
\value{
A character vector of Internet Archive item IDs.
}
\description{
Perform an advanced search of the Internet Archive, specifying which metadata
fields to search. Note that all searches are in the form of "contains," i.e.,
the title contains the search term.
}
\examples{
query1 <- c("title" = "damnation of theron ware")
ia_search(query1)
query2 <- c("title" = "damnation of theron ware",
            "contributor" = "gutenberg")
ia_search(query2)
}
\references{
See the documentation on the Internet Archive's
  \href{https://archive.org/advancedsearch.php}{advanced search page}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ia_keyword_search.R
\name{ia_keyword_search}
\alias{ia_keyword_search}
\title{Perform an simple keyword search of the Internet Archive.}
\usage{
ia_keyword_search(keywords, num_results = 5, page = 1, print_total = TRUE)
}
\arguments{
\item{keywords}{The keywords to search for.}

\item{num_results}{The number of results to return per page.}

\item{page}{When results are paged, which page of results to return.}

\item{print_total}{Should the total number of results for this query be
printed as a message?}
}
\value{
A character vector of Internet Archive item IDs.
}
\description{
Perform an simple keyword search of the Internet Archive.
}
\examples{
ia_keyword_search("isaac hecker", num_results = 20)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ia_download.R
\name{ia_download}
\alias{ia_download}
\title{Download files for Internet Archive items.}
\usage{
ia_download(files, dir = ".", extended_name = TRUE, overwrite = FALSE,
  silence = FALSE)
}
\arguments{
\item{files}{A data frame of files returned by \link{ia_files}. You should
filter this data frame to download only the files that you actually want.}

\item{dir}{The directory in which to save the downloaded files.}

\item{extended_name}{If this argument is \code{FALSE}, then the downloaded
file will have a filename in the following format:
\code{itemidentifier.extension}, e.g., \code{thedamnationofth00133gut.txt}.
If there are multiple files of the same file type for an item, then the
file names will not be unique. If this argument is \code{TRUE}, them the
downloaded file will have a filename in the following format:
\code{itemidentifier-original-filename.extension}, e.g.,
\code{thedamnationofth00133gut-133.txt}.}

\item{overwrite}{If \code{TRUE}, this function will download all files and
overwrite them on disk if they have already been downloaded. If
\code{FALSE}, then if a file already exists on disk it will not be
downloaded again but other downloads will proceed normally.}

\item{silence}{If false, print the item IDs as they are downloaded.}
}
\value{
A data frame including the file names of the downloaded files.
}
\description{
Download files for Internet Archive items.
}
\examples{
\dontrun{
if (require(dplyr)) {
  dir <- tempdir()
  ia_get_items("thedamnationofth00133gut") \%>\%
    ia_files() \%>\%
    filter(type == "txt") \%>\% # download only the files we want
    ia_download(dir = dir, extended_name = FALSE)
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ia_list_fields.R
\name{ia_list_fields}
\alias{ia_list_fields}
\title{List accepted metadata fields}
\usage{
ia_list_fields()
}
\value{
A list of the accepted metadata fields
}
\description{
List accepted metadata fields
}
\examples{
ia_list_fields()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ia_get_item.R
\name{ia_get_items}
\alias{ia_get_items}
\title{Get the metadata for Internet Archive items}
\usage{
ia_get_items(item_id, silence = FALSE)
}
\arguments{
\item{item_id}{A character vector containing the ID for an Internet Archive
item. This argument is vectorized, so you can retrieve multiple items at
once.}

\item{silence}{If false, print the item IDs as they are retrieved.}
}
\value{
A list containing the metadata returned by the API. List names
  correspond to the item IDs.
}
\description{
Get the metadata for Internet Archive items
}
\examples{
\dontrun{
ia_get_items("thedamnationofth00133gut")

ats_query <- c("publisher" = "american tract society")
ids       <- ia_search(ats_query, num_results = 2)
ia_get_items(ids)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ia_browse.R
\name{ia_browse}
\alias{ia_browse}
\title{Open an Internet Archive item in the browser}
\usage{
ia_browse(item_id, type = c("details", "stream"))
}
\arguments{
\item{item_id}{The item identifier. If multiple item identifiers are passed
in, only the first will be opened.}

\item{type}{Which page to open: \code{details} is the metadata page,
\code{stream} is the viewing page for items which are associated with a
PDF.}
}
\value{
Returns the item ID(s) passed to the function.
}
\description{
Open an Internet Archive item in the browser
}
\examples{
# Distinguished Converts to Rome in America
ia_browse("distinguishedcon00scanuoft")
}
