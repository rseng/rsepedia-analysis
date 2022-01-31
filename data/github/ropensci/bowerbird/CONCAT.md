
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

[![R-CMD-check](https://github.com/ropensci/bowerbird/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/bowerbird/actions)
[![Codecov test
coverage](https://codecov.io/gh/ropensci/bowerbird/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/bowerbird?branch=master)
[![](https://badges.ropensci.org/139_status.svg)](https://github.com/ropensci/onboarding/issues/139)
<!-- badges: end -->

# Bowerbird

<img align="right" src="https://rawgit.com/ropensci/bowerbird/master/inst/extdata/bowerbird.svg" />

Often it’s desirable to have local copies of third-party data sets.
Fetching data on the fly from remote sources can be a great strategy,
but for speed or other reasons it may be better to have local copies.
This is particularly common in environmental and other sciences that
deal with large data sets (e.g. satellite or global climate model
products). Bowerbird is an R package for maintaining a local collection
of data sets from a range of data providers.

A comprehensive introduction to bowerbird can be found at
<https://docs.ropensci.org/bowerbird/articles/bowerbird.html>, along
with full package documentation.

## Installing

``` r
## use the SCAR r-universe package repository
options(repos = c(SCAR = "https://scar.r-universe.dev", CRAN = "https://cloud.r-project.org"))

## install
install.packages("bowerbird")

## or install from github
##install.packages("remotes") ## if needed
remotes::install_github("ropensci/bowerbird", build_vignettes = TRUE)
```

## Usage overview

### Configuration

Build up a configuration by first defining global options such as the
destination on your local file system. Commonly you would choose this
destination data directory to be a persistent location, suitable for a
data library. For demonstration purposes here we’ll just use a temporary
directory::

``` r
library(bowerbird)
my_directory <- tempdir()
cf <- bb_config(local_file_root = my_directory)
```

Bowerbird must then be told which data sources to synchronize. Let’s use
data from the Australian 2016 federal election, which is provided as one
of the example data sources:

``` r
my_source <- bb_example_sources("Australian Election 2016 House of Representatives data")

## add this data source to the configuration
cf <- bb_add(cf, my_source)
```

Once the configuration has been defined and the data source added to it,
we can run the sync process. We set `verbose = TRUE` here so that we see
additional progress output:

``` r
status <- bb_sync(cf, verbose = TRUE)
```

    ##  
    ## Tue Nov 30 16:32:35 2021 
    ## Synchronizing dataset: Australian Election 2016 House of Representatives data 
    ## Source URL http://results.aec.gov.au/20499/Website/HouseDownloadsMenu-20499-Csv.htm 
    ## -------------------------------------------------------------------------------------------- 
    ##  
    ##  this dataset path is: /tmp/data/results.aec.gov.au/20499/Website 
    ##  visiting http://results.aec.gov.au/20499/Website/HouseDownloadsMenu-20499-Csv.htm ... done. 
    ##  downloading file 1 of 47: http://results.aec.gov.au/20499/Website/Downloads/HouseCandidatesDownload-20499.csv ...  done. 
    ##  downloading file 2 of 47: http://results.aec.gov.au/20499/Website/Downloads/HouseMembersElectedDownload-20499.csv ...  done. 
    ##  downloading file 3 of 47: http://results.aec.gov.au/20499/Website/Downloads/HouseNominationsByStateDownload-20499.csv ...  done. 
    ##  
    ##  [... output truncated] 
    ##  
    ## Tue Nov 30 16:32:49 2021 dataset synchronization complete: Australian Election 2016 House of Representatives data

Congratulations\! You now have your own local copy of your chosen data
set. This particular example is fairly small (about 10MB), so it should
not take too long to download. Details of the files in this data set are
given in the `status$files` object:

``` r
status$files
## [[1]]
## # A tibble: 47 × 3
##    url                                file                               note   
##    <chr>                              <chr>                              <chr>  
##  1 http://results.aec.gov.au/20499/W… /tmp/data/results.aec.gov.au/2049… downlo…
##  2 http://results.aec.gov.au/20499/W… /tmp/data/results.aec.gov.au/2049… downlo…
##  3 http://results.aec.gov.au/20499/W… /tmp/data/results.aec.gov.au/2049… downlo…
##  4 http://results.aec.gov.au/20499/W… /tmp/data/results.aec.gov.au/2049… downlo…
##  5 http://results.aec.gov.au/20499/W… /tmp/data/results.aec.gov.au/2049… downlo…
##  6 http://results.aec.gov.au/20499/W… /tmp/data/results.aec.gov.au/2049… downlo…
##  7 http://results.aec.gov.au/20499/W… /tmp/data/results.aec.gov.au/2049… downlo…
##  8 http://results.aec.gov.au/20499/W… /tmp/data/results.aec.gov.au/2049… downlo…
##  9 http://results.aec.gov.au/20499/W… /tmp/data/results.aec.gov.au/2049… downlo…
## 10 http://results.aec.gov.au/20499/W… /tmp/data/results.aec.gov.au/2049… downlo…
## # … with 37 more rows
```

At a later time you can re-run this synchronization process. If the
remote files have not changed, and assuming that your configuration has
the `clobber` parameter set to 0 (“do not overwrite existing files”) or
1 (“overwrite only if the remote file is newer than the local copy”)
then the sync process will run more quickly because it will not need to
re-download any data files.

## Data source definitions

The [blueant](https://github.com/AustralianAntarcticDivision/blueant)
package provides a suite of bowerbird data source definitions themed
around Southern Ocean and Antarctic data, and includes a range of
oceanographic, meteorological, topographic, and other environmental data
sets.

## Other packages

Many other data-retrieval R packages exist. bowerbird is perhaps most
similar to the
[rdataretriever](https://cran.r-project.org/package=rdataretriever).
This package provides an R interface to the (Python-based) [Data
Retriever](http://www.data-retriever.org/), which in turn provides (at
time of writing) access to 85 ecological data sets. A quick comparison:

### rdataretriever

  - requires `retriever` to be installed, either as a Python package or
    via a platform-specific installer (see
    <http://www.data-retriever.org/>)

  - makes efforts to clean and standardize the data that it downloads,
    and get them into a consistent format on the user’s system

  - designed to make it easy for users to get on with the business of
    using those data sets

  - carries the tradeoff that adding new data sets (and maintaining the
    existing ones) takes a bit of effort, and it can be cumbersome to
    deal with data sets that contain many files, particularly if new
    files get added on a regular basis (e.g. satellite environmental
    data).

### bowerbird

  - pure R, no other system dependencies

  - designed to make it easy for users to keep a local, up-to-date
    collection of files from remote providers. It can do recursive
    downloads, and so is particularly suitable for collections that are
    structured as a large number of individual files in yearly or other
    subdirectories (typical of e.g. satellite or climate model data)

  - simply mirrors remote data to your local system, without attempting
    to reformat the data files or do anything else clever with them
    (other than uncompress, if needed). It just grabs them and saves
    them in whatever format the provider uses

  - the upside is that it is intended to be easy to write bowerbird
    definitions for new data sources. In many cases, it is only
    necessary to specify some metadata and the top-level URL, and
    bowerbird can recursively download linked resources from there

  - bowerbird itself contains only a few example data sets, but data
    definitions are available from other packages
    (e.g. [blueant](https://github.com/AustralianAntarcticDivision/blueant),
    \~55 marine/Southern Ocean data sets).

The rdataretriever and bowerbird packages are both part of the rOpenSci
project.

[![ropensci\_footer](https://ropensci.org/public_images/scar_footer.png)](https://ropensci.org)
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
# Contributions

Suggestions, bug reports, and code are welcome. Please [fork this repository](https://help.github.com/articles/fork-a-repo/) and submit a pull request.

Requests are also welcome: open an [issue](https://github.com/ropensci/bowerbird/issues).

# Code of Conduct

Maintainers and contributors must follow this repository's [code of conduct](CODE_OF_CONDUCT.md).

---
output: github_document
editor_options: 
  chunk_output_type: console
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  ##comment = "#>",
  fig.path = "README-"
)
options(tibble.width = 80, tibble.print_max = 10, tibble.print_min = 10)
```

<!-- badges: start -->
[![R-CMD-check](https://github.com/ropensci/bowerbird/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/bowerbird/actions)
[![Codecov test coverage](https://codecov.io/gh/ropensci/bowerbird/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/bowerbird?branch=master)
[![](https://badges.ropensci.org/139_status.svg)](https://github.com/ropensci/onboarding/issues/139)
<!-- badges: end -->

# Bowerbird

<img align="right" src="https://rawgit.com/ropensci/bowerbird/master/inst/extdata/bowerbird.svg" />

Often it's desirable to have local copies of third-party data sets. Fetching data on the fly from remote sources can be a great strategy, but for speed or other reasons it may be better to have local copies. This is particularly common in environmental and other sciences that deal with large data sets (e.g. satellite or global climate model products). Bowerbird is an R package for maintaining a local collection of data sets from a range of data providers.

A comprehensive introduction to bowerbird can be found at https://docs.ropensci.org/bowerbird/articles/bowerbird.html, along with full package documentation.

## Installing

```{r eval = FALSE}
## use the SCAR r-universe package repository
options(repos = c(SCAR = "https://scar.r-universe.dev", CRAN = "https://cloud.r-project.org"))

## install
install.packages("bowerbird")

## or install from github
##install.packages("remotes") ## if needed
remotes::install_github("ropensci/bowerbird", build_vignettes = TRUE)

```

## Usage overview

### Configuration

Build up a configuration by first defining global options such as the destination on your local file system. Commonly you would choose this destination data directory to be a persistent location, suitable for a data library. For demonstration purposes here we'll just use a temporary directory::

```{r include = FALSE}
## code not shown in the README: (1) use load_all() instead of library() for convenience, (2) use anonymous "temporary" dir, and (3) hide the progress bar
devtools::load_all()
my_data_dir <- "/tmp/data"
if (!dir.exists(my_data_dir)) dir.create(my_data_dir)
if (dir.exists(file.path(my_data_dir, "results.aec.gov.au"))) unlink(file.path(my_data_dir, "results.aec.gov.au"), recursive = TRUE)
cf <- bb_config(local_file_root = my_data_dir)
mysrc <- bb_example_sources("Australian Election 2016 House of Representatives data") %>% bb_modify_source(method = list(show_progress = FALSE))
cf <- cf %>% bb_add(mysrc)
```


```{r eval = FALSE}
library(bowerbird)
my_directory <- tempdir()
cf <- bb_config(local_file_root = my_directory)
```

Bowerbird must then be told which data sources to synchronize. Let's use data from the Australian 2016 federal election, which is provided as one of the example data sources:

```{r eval = FALSE}
my_source <- bb_example_sources("Australian Election 2016 House of Representatives data")

## add this data source to the configuration
cf <- bb_add(cf, my_source)
```

Once the configuration has been defined and the data source added to it, we can run the sync process. We set `verbose = TRUE` here so that we see additional progress output:

```{r eval = FALSE}
status <- bb_sync(cf, verbose = TRUE)
```

```{r echo = FALSE, message = FALSE}
## code not shown in README: capture the output and trim it down a bit
op <- capture.output(status <- bb_sync(cf, verbose = TRUE))
idx <- grepl("^ downloading file", op)
if (sum(idx)>5) {
  op[which(idx)[4]] <- ""
  op[which(idx)[5]] <- " [... output truncated]"
  idx[which(idx)[1:5]] <- FALSE
  op <- op[!idx]
  }
#op ## gives ## [1] "" etc
for (oo in op) cat(oo, "\n")
```

Congratulations! You now have your own local copy of your chosen data set. This particular example is fairly small (about 10MB), so it should not take too long to download. Details of the files in this data set are given in the `status$files` object:

```{r}
status$files
```

At a later time you can re-run this synchronization process. If the remote files have not changed, and assuming that your configuration has the `clobber` parameter set to 0 ("do not overwrite existing files") or 1 ("overwrite only if the remote file is newer than the local copy") then the sync process will run more quickly because it will not need to re-download any data files.


## Data source definitions

The [blueant](https://github.com/AustralianAntarcticDivision/blueant) package provides a suite of bowerbird data source definitions themed around Southern Ocean and Antarctic data, and includes a range of oceanographic, meteorological, topographic, and other environmental data sets.


## Other packages

Many other data-retrieval R packages exist. bowerbird is perhaps most similar to the [rdataretriever](https://cran.r-project.org/package=rdataretriever). This package provides an R interface to the (Python-based) [Data Retriever](http://www.data-retriever.org/), which in turn provides (at time of writing) access to 85 ecological data sets. A quick comparison:

### rdataretriever

- requires `retriever` to be installed, either as a Python package or via a platform-specific installer (see http://www.data-retriever.org/)

- makes efforts to clean and standardize the data that it downloads, and get them into a consistent format on the user's system

- designed to make it easy for users to get on with the business of using those data sets

- carries the tradeoff that adding new data sets (and maintaining the existing ones) takes a bit of effort, and it can be cumbersome to deal with data sets that contain many files, particularly if new files get added on a regular basis (e.g. satellite environmental data).

### bowerbird

- pure R, no other system dependencies

- designed to make it easy for users to keep a local, up-to-date collection of files from remote providers. It can do recursive downloads, and so is particularly suitable for collections that are structured as a large number of individual files in yearly or other subdirectories (typical of e.g. satellite or climate model data)

- simply mirrors remote data to your local system, without attempting to reformat the data files or do anything else clever with them (other than uncompress, if needed). It just grabs them and saves them in whatever format the provider uses

- the upside is that it is intended to be easy to write bowerbird definitions for new data sources. In many cases, it is only necessary to specify some metadata and the top-level URL, and bowerbird can recursively download linked resources from there

- bowerbird itself contains only a few example data sets, but data definitions are available from other packages (e.g. [blueant](https://github.com/AustralianAntarcticDivision/blueant), ~55 marine/Southern Ocean data sets).

The rdataretriever and bowerbird packages are both part of the rOpenSci project.

[![ropensci_footer](https://ropensci.org/public_images/scar_footer.png)](https://ropensci.org)

---
title: "Data provenance"
author: "Ben Raymond, Michael Sumner"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Data provenance}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

Bowerbird will maintain a local collection of data files, sourced from external data providers. When one does an analysis using such files, it is useful to know *which* files were used, and the *provenance* of those files. This information will assist in making analyses reproducible.

Bowerbird itself actually knows very little about the data files that it maintains, particularly if it is using the `bb_handler_rget` method. In this case, for a given data source it typically knows the URL and flags to pass to `bb_rget`, and some basic metadata (primarily intended to be read by the user). This is by design: the heavy lifting involved in mirroring a remote data source is generally handballed `bb_rget`. This simplifies both the bowerbird code as well as the process of writing and maintaining data sources.

The downside of this is that bowerbird can't necessarily answer questions about data provenance in specific detail. Say that we are using data from a data source that contains many separate files, one per day, spanning many years, and we run an analysis that uses only files from the year 1999. Ideally we would like a function that could tell us exactly which files from the complete data set are needed for our particular analysis, but this is impossible without specific knowledge of how the data source is structured. Analogous situations exist with data sources that split geographic space across files, or split different parameters across files.

Despite this, bowerbird can provide some general information to assist with data provenance issues.

## Which files were used in an analysis

Say we have a bowerbird config:
```{r}
library(bowerbird)
cf <- bb_config("/some/local/path") %>% bb_add(bb_example_sources())
```

We can ask bowerbird where these files are stored locally:
```{r}
source_dirs <- bb_data_source_dir(cf)
knitr::kable(data.frame(name = cf$data_sources$name, source_dir = source_dirs))
```

The directory for each data source will contain *all* of the files associated with that data source, not just the files used in a particular analysis. (In extreme cases this directory might even contain files from other data sources, although this should only happen if multiple data sources have overlapping `source_url` paths, which ought not to be a common occurrence.)

So while this directory is not the *minimal* set of files needed to reproduce an analysis, it does at least contain that set, and could be used to store those files in an online repository that assigns DOIs (such as figshare, see https://cran.r-project.org/package=rfigshare, or zenodo), or bundle into a docker image (see e.g https://github.com/o2r-project/containerit).

The complete list of files within a data source's directory can be generated using a standard directory listing (e.g. `list.files(path = my_source_dir, recursive = TRUE`), but see also the `bb_fingerprint` function described below, which provides additional information about the files.

### Refining the list of files

Ascertaining exactly which files were used in a particular analysis is a task that is better handled by the code being used to do the file-reading and analysis, not by the repository-management (bowerbird) code.

We are unaware of any general solutions that will keep track of the files used by an arbitrary chunk of R code.

The recordr package (https://github.com/NCEAS/recordr) comes close: it will collect info about which files were read or written. However, this only works for file types that have read functions implemented in recordr. At the time of writing, this does not cover netcdf files or other files typically used for environmental data, and so its application in that domain is likely to be of limited value.

Packages that provide specific data access functionality (e.g. https://github.com/AustralianAntarcticDivision/raadtools) might provide mechanisms for tracking which data files are needed in order to fulfil particular data queries.

## The provenance of files used in an analysis

Provenance: where the files came from, when they were downloaded, their version, etc.

The `bb_fingerprint` function, given a data repository configuration, will return the timestamp of download and hashes of all files associated with its data sources. Thus, for all of these files, we have established where they came from (the data source ID), when they were downloaded, and a hash so that later versions of those files can be compared to detect changes.

### A word on digital object identifiers

When a data source is defined, it should have used its DOI as its data source ID (if it has a DOI, otherwise some other unique identifier). The idea of a DOI is that a data set that has a DOI should be peristent (accessible for long term use), and the DOI should uniquely identify the data resource that it is assigned to. When a data set changes in a substantial manner, and/or it is necessary to identify both the original and the changed material, a new DOI should be assigned. Thus, knowing the DOI gives some indication of the data provenance.

However, it is worth noting that the DOI might not uniquely identify the state of the data source. For example, the NSIDC SMMR-SSM/I Nasateam near-real-time sea ice concentration data set (DOI http://doi.org/10.5067/U8C09DWVX9LM) is updated each day (new files are added) as new data is acquired by the satellites. At periodic intervals, these files are subjected to more rigorous quality control and post-processing, and then moved into a different data set (http://doi.org/10.5067/8GQ8LZQVL0VL). Thus, knowing the DOI of the near-real-time data set does not uniquely identify the files that were included in that data set at a given point in time.

Similar ambiguity can arise for other reasons, including data corrections. Say that one data file within a large collection was found to have errors and was corrected: it is at the discretion of the data provider whether the data set as a whole receives a new DOI because of this change.

---
title: "Bowerbird"
author: "Ben Raymond, Michael Sumner"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
    number_sections: true
vignette: >
  %\VignetteIndexEntry{bowerbird}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

# Bowerbird

### Table of contents

<ul>
<li><a href="#overview">Overview</a></li>
<li><a href="#installing">Installing</a></li>
<li><a href="#usage-overview">Usage overview</a></li>
<li><a href="#users-level-of-usage-and-expected-knowledge">Users: level of usage and expected knowledge</a></li>
<li><a href="#defining-data-sources">Defining data sources</a></li>
<li><a href="#nuances">Nuances</a></li>
<li><a href="#developer-notes">Developer notes</a></li>
<li><a href="#data-source-summary">Data source summary</a></li>
</ul>
&nbsp;

## Overview

Often it's desirable to have local copies of third-party data sets. Fetching data on the fly from remote sources can be a great strategy, but for speed or other reasons it may be better to have local copies. This is particularly common in environmental and other sciences that deal with large data sets (e.g. satellite or global climate model products). Bowerbird is an R package for maintaining a local collection of data sets from a range of data providers.

Bowerbird can be used in several different modes:

- interactively from the R console, to download or update data files on an as-needed basis
- from the command line, perhaps as a regular scheduled task
- programatically, including from within other R packages, scripts, or R markdown documents that require local copies of particular data files.

When might you consider using bowerbird rather than, say, [curl](https://cran.r-project.org/package=curl) or [crul](https://cran.r-project.org/package=crul)? The principal advantage of bowerbird is that it can download files recursively. In many cases, it is only necessary to specify the top-level URL, and bowerbird can recursively download linked resources. Bowerbird can also:


- decompress downloaded files (if the remote server provides them in, say, zipped or gzipped form).

- incrementally update files that you have previously downloaded. Bowerbird can be instructed not to re-download files that exist locally, unless they have changed on the remote server. Compressed files will also only be decompressed if changed.


## Installing

```{r eval = FALSE}
install.packages("remotes")
remotes::install_github("ropensci/bowerbird", build_vignettes = TRUE)
```

## Usage overview

### Configuration

Build up a configuration by first defining global options such as the destination on your local file system:

```{r eval = FALSE}
library(bowerbird)
my_directory <- "~/my/data/directory"
cf <- bb_config(local_file_root = my_directory)
```

Bowerbird must then be told which data sources to synchronize. Let's use data from the Australian 2016 federal election, which is provided as one of the example data sources:

```{r eval = FALSE}
my_source <- bb_example_sources("Australian Election 2016 House of Representatives data")

## add this data source to the configuration
cf <- bb_add(cf, my_source)
```

Once the configuration has been defined and the data source added to it, we can run the sync process. We set `verbose = TRUE` here so that we see additional progress output:

```{r eval = FALSE}
status <- bb_sync(cf, verbose = TRUE)
```

Congratulations! You now have your own local copy of your chosen data set. This particular example is fairly small (about 10MB), so it should not take too long to download. Details of the files in this data set are given in the `status$files` object:

```{r eval = FALSE}
status$files
```

At a later time you can re-run this synchronization process. If the remote files have not changed, and assuming that your configuration has the `clobber` parameter set to 0 ("do not overwrite existing files") or 1 ("overwrite only if the remote file is newer than the local copy") then the sync process will run more quickly because it will not need to re-download any data files.


## Users: level of usage and expected knowledge

Users can interact with bowerbird at several levels, with increasing levels of complexity:

1. **Using bowerbird with data source definitions that have been written by someone else**. This is fairly straightforward. 

1. **Writing your own data source definitions so that you can download data from a new data provider, but using an existing handler such as `bb_handler_rget`**. This is a little more complicated. You will need reasonable knowledge of how your data provider disseminates its files (including e.g. the source URL from which data files are to be downloaded, and how the data repository is structured). Be prepared to fiddle with `rget` settings to accommodate provider-specific requirements (e.g. controlling recursion behaviour). The "Defining data sources" section below provides guidance and examples on writing data source definitions.

1. **Writing your own handler function for data providers that do not play nicely with the packaged handlers (`bb_handler_rget`, `bb_handler_oceandata`, `bb_handler_earthdata`)**. This is the trickiest, and at the time of writing we have not provided much guidance on how to do this. See the "Writing new data source handlers" section, below.

It is expected that the majority of users will fall into one of the first two categories.

## Defining data sources

### Prepackaged data source definitions

A few example data source definitions are provided as part of the bowerbird package --- see `bb_example_sources()` (these are also listed at the bottom of this document).

Other packages provide data source definitions that can be used with bowerbird. The [blueant](https://github.com/AustralianAntarcticDivision/blueant) package provides a suite of bowerbird data source definitions themed around Southern Ocean and Antarctic data, and includes a range of oceanographic, meteorological, topographic, and other environmental data sets.

### Defining your own data sources

The general bowerbird workflow is to build up a configuration with one or more data sources, and pass that configuration object to the `bb_sync` function to kick off the download process. Each data source contains the details required by bowerbird to be able to fetch it, including a *handler* function that bb_sync will call to do the actual download.

The `bb_handler_rget` function is a generic handler function that will be suitable for many data sources. Note that this `bb_handler_rget` function is not intended to be called directly by the user, but is specified as part of a data source specification. The `bb_sync` function calls `bb_handler_rget` during its run, passing appropriate parameters as it does so.

`bb_handler_rget` is a wrapper around `bb_rget`, which is a recursive file downloading utility. Typically, one only needs to define a data source in terms of its top-level URL and appropriate flags to pass to `bb_rget`, along with some basic metadata (primarily intended to be read by the user).

Specifying a data source is done by the `bb_source` function. This can seem a little daunting, so let's work through some examples. Most of these examples are included in `bb_example_sources()`.

#### Example 1: a single data file

Say we've found [this bathymetric data set](https://doi.org/10.4225/25/53D9B12E0F96E) and we want to define a data source for it. It's a single zip file that contains some ArcGIS binary grids. A minimal data source definition might look like this:

```{r eval = FALSE}
src1 <- bb_source(
    name = "Geoscience Australia multibeam bathymetric grids of the Macquarie Ridge",
    id = "10.4225/25/53D9B12E0F96E",
    doc_url = "https://doi.org/10.4225/25/53D9B12E0F96E",
    license = "CC-BY 4.0",
    citation = "Spinoccia, M., 2012. XYZ multibeam bathymetric grids of the Macquarie Ridge. Geoscience Australia, Canberra.",
    source_url = "http://www.ga.gov.au/corporate_data/73697/Macquarie_ESRI_Raster.zip",
    method = list("bb_handler_rget"))
```

The parameters provided here are all mandatory:


- `id` is a unique identifier for the dataset, and should be something that changes when the data set is updated. Its DOI, if it has one, is ideal for this. Otherwise, if the original data provider has an identifier for this dataset, that is probably a good choice here (include the data version number if there is one)
- `name` is a human-readable but still unique identifier
- `doc_url` is a link to a metadata record or documentation page that describes the data in detail
- `license` is the license under which the data are being distributed, and is required so that users are aware of the conditions that govern the usage of the data
- `citation` gives citation details for the data source. It's generally considered good practice to cite data providers, and indeed under some data licenses this is in fact mandatory
- the `method` parameter is specified as a list, where the first entry is the name of the handler function that will be used to retrieve this data set (`bb_handler_rget`, in this case) and the remaining entries are data-source-specific arguments to pass to that function (none here)
- and finally the `source_url` tells bowerbird where to go to get the data.


Add this data source to a configuration and synchronize it:
```{r eval = FALSE}
cf <- bb_config("c:/temp/data/bbtest") %>% bb_add(src1)
status <- bb_sync(cf)
```

This should have caused the zip file to be downloaded to your local machine, in this case into the `c:/temp/data/bbtest/www.ga.gov.au/corporate_data/73697` directory. 

There are a few additional entries that we might consider for this data source, particularly if we are going to make it available for other users.

Firstly, having the zip file locally is great, but we will need to unzip it before we can actually use it. Bowerbird can do this by specifying a `postprocess` action:

```{r eval = FALSE}
src1 <- bb_source(..., postprocess = list("bb_unzip"))
```

For the benefit of other users, we might also add the `description`, `collection_size`, and `data_group` parameters:

- `description` provides a plain-language description of the data set, so that users can get an idea of what it contains (for full details they can consult the `doc_url` link that you already provided)
- `collection_size` is the approximate disk space (in GB) used by the data collection. Some collections are very large! This parameter obviously gives an indication of the storage space required, but also the download size (noting though that some data sources deliver compressed files, so the download size might be much smaller)
- `data_group` is a descriptive or thematic group name that this data set belongs to. This can also help users find data sources of interest to them
- `access_function` can be used to suggest to users an appropriate function to read these data files. In this case the files can be read by the `raster` function (from the `raster` package).

So our full data source definition now looks like:

```{r eval = FALSE}
src1 <- bb_source(
    name = "Geoscience Australia multibeam bathymetric grids of the Macquarie Ridge",
    id = "10.4225/25/53D9B12E0F96E",
    description = "This is a compilation of all the processed multibeam bathymetry data that are publicly available in Geoscience Australia's data holding for the Macquarie Ridge.",
    doc_url = "https://doi.org/10.4225/25/53D9B12E0F96E",
    license = "CC-BY 4.0",
    citation = "Spinoccia, M., 2012. XYZ multibeam bathymetric grids of the Macquarie Ridge. Geoscience Australia, Canberra.",
    source_url = "http://www.ga.gov.au/corporate_data/73697/Macquarie_ESRI_Raster.zip",
    method = list("bb_handler_rget"),
    postprocess = list("bb_unzip"),
    collection_size = 0.4,
    access_function = "raster::raster",
    data_group = "Topography")
```

#### Example 2: multiple files via recursive download

This data source (Australian Election 2016 House of Representatives data) is provided as one of the example data sources in `bb_example_sources()`, but let's look in a little more detail here.

The primary entry point to this data set is an HTML index page, which links to a number of data files. In principle we could generate a list of all of these data files and download them one by one, but that would be tedious and prone to breaking (if the data files changed our hard-coded list would no longer be correct). Instead we can start at the HTML index and recursively download linked data files.

The definition for this data source is:

```{r eval = FALSE}
src2 <- bb_source(
    name = "Australian Election 2016 House of Representatives data",
    id = "aus-election-house-2016",
    description = "House of Representatives results from the 2016 Australian election.",
    doc_url = "http://results.aec.gov.au/",
    citation = "Copyright Commonwealth of Australia 2017. As far as practicable, material for which the copyright is owned by a third party will be clearly labelled. The AEC has made all reasonable efforts to ensure that this material has been reproduced on this website with the full consent of the copyright owners.",
    source_url = c("http://results.aec.gov.au/20499/Website/HouseDownloadsMenu-20499-Csv.htm"),
    license = "CC-BY",
    method = list("bb_handler_rget", level = 1, accept_download = "\\.csv$"),
    collection_size=0.01)
```

Most of these parameters will be familiar from the previous example, but the `method` definition is more complex. Let's look at the entries in the `method` list (these are all parameters that get passed to the `bb_rget()` function):

- `level = 1` specifies that `bb_rget` should download recursively, but only recurse down one level (i.e. follow links found in the top-level `source_url` document, but don't recurse any deeper. If, say, we specified `level = 2`, then `bb_rget` would follow links from the top-level document as well as links found in those linked documents.)
- `accept_download = "\\.csv$"` means that we only want to retrieve csv files.

Add this data source to a configuration and synchronize it:
```{r eval = FALSE}
cf <- bb_config("c:/temp/data/bbtest") %>% bb_add(src2)
status <- bb_sync(cf)
```

Once again the data have been saved into a subdirectory that reflects the URL structure (`c:/temp/data/bbtest/results.aec.gov.au/20499/Website/Downloads`). If you examine that directory, you will see that it contains around 50 separate csv files, each containing a different component of the data set.

You can immediately see that by using a recursive download, not only did we not need to individually specify all 50 of those data files, but if the data provider adds new files in the future the recursive download process will automatically find them (so long as they are linked from the top-level `source_url` document).

#### Example 3: an Earthdata source

The [Earthdata system](https://earthdata.nasa.gov/) is one of NASA's data management systems and home to a vast range of Earth science data from satellites, aircraft, field measurements, and other sources. Say you had a rummage through their [data catalogue](https://search.earthdata.nasa.gov/) and found yourself wanting a copy of [Sea Ice Trends and Climatologies from SMMR and SSM/I-SSMIS](https://doi.org/10.5067/IJ0T7HFHB9Y6).

Data sources served through the Earthdata system require users to have an Earthdata account, and to log in with their credential when downloading data. Bowerbird's `bb_handler_earthdata` function eases some of the hassle involved with these Earthdata sources.

First, let's create an account and get ourselves access to the data.

1. create an Earthdata login via https://wiki.earthdata.nasa.gov/display/EL/How+To+Register+With+Earthdata+Login if you don't already have one

1. we need to know the URL of the actual data. The [metadata record](https://doi.org/10.5067/IJ0T7HFHB9Y6) for this data set contains a "Download via HTTPS" link, which in turn directs the user to this URL: https://daacdata.apps.nsidc.org/pub/DATASETS/nsidc0192_seaice_trends_climo_v3/. That's the data URL (which will be used as the `source_url` in our data source definition)

1. browse to the [that data URL](https://daacdata.apps.nsidc.org/pub/DATASETS/nsidc0192_seaice_trends_climo_v3/), using your Earthdata login to authenticate. When you use access an Earthdata application for the first time, you will be requested to authorize it so that it can access data using your credentials (see https://wiki.earthdata.nasa.gov/display/EL/How+To+Register+With+Earthdata+Login). This dataset is served by the NSIDC DAAC application, so you will need to authorize this application (either through browsing as just described, or go to 'My Applications' at https://urs.earthdata.nasa.gov/profile and add the application 'nsidc-daacdata' to your list of authorized applications)

You only need to create an Earthdata login once. If you want to download other Earthdata data sets via bowerbird, you'll use the same credentials, but note that you may need to authorize additional applications, depending on where your extra data sets come from.

Now that we have access to the data, we can write our bowerbird data source:

```{r eval = FALSE}
src3 <- bb_source(
    name = "Sea Ice Trends and Climatologies from SMMR and SSM/I-SSMIS, Version 3",
    id = "10.5067/IJ0T7HFHB9Y6",
    description = "NSIDC provides this data set to aid in the investigations of the variability and trends of sea ice cover. Ice cover in these data are indicated by sea ice concentration: the percentage of the ocean surface covered by ice. The ice-covered area indicates how much ice is present; it is the total area of a pixel multiplied by the ice concentration in that pixel. Ice persistence is the percentage of months over the data set time period that ice existed at a location. The ice-extent indicates whether ice is present; here, ice is considered to exist in a pixel if the sea ice concentration exceeds 15 percent. This data set provides users with data about total ice-covered areas, sea ice extent, ice persistence, and monthly climatologies of sea ice concentrations.",
    doc_url = "https://doi.org/10.5067/IJ0T7HFHB9Y6",
    citation = "Stroeve J, Meier WN (2018) Sea Ice Trends and Climatologies from SMMR and SSM/I-SSMIS, Version 3. [Indicate subset used]. Boulder, Colorado USA. NASA National Snow and Ice Data Center Distributed Active Archive Center. doi:10.5067/IJ0T7HFHB9Y6. [Date Accessed].",
    source_url = c("https://daacdata.apps.nsidc.org/pub/DATASETS/nsidc0192_seaice_trends_climo_v3/"),
    license = "Please cite, see http://nsidc.org/about/use_copyright.html",
    authentication_note = "Requires Earthdata login, see https://wiki.earthdata.nasa.gov/display/EL/How+To+Register+With+Earthdata+Login . Note that you will also need to authorize the application 'nsidc-daacdata' (see 'My Applications' at https://urs.earthdata.nasa.gov/profile)",
    method = list("bb_handler_earthdata", level = 4, relative = TRUE, accept_download = "\\.(s|n|png|txt)$"),
    user = "your_earthdata_username",
    password = "your_earthdata_password",
    collection_size = 0.02,
    data_group = "Sea ice")
```

This is very similar to our previous examples, with these differences:

- the `method` specifies `bb_handler_earthdata` (whereas previously we used `bb_handler_rget`). The `bb_handler_earthdata` is actually just a wrapper around `bb_handler_rget`, but it takes care of some Earthdata-specific details like authentication using your Earthdata credentials
- `relative = TRUE` means that `bb_rget` will only follow relative links (i.e. links of the form `<a href="/some/directory/">...</a>`, which by definition must be on the same server as our `source_url`). Absolute links (i.e. links of the form `<a href="http://some.other.server/some/path">...</a>` will not be followed. This is another mechanism to prevent the recursive download from downloading stuff we don't want
- we have specified the file types that we want via the `accept_download` parameter.


Note that if you were providing this data source definition for other people to use, you would obviously not want to hard-code your Earthdata credentials in the `user` and `password` slots. In this case, specify the credentials as empty strings, and also include `warn_empty_auth = FALSE` in the data source definition (this suppresses the warning that `bb_source` would otherwise give you about missing credentials):

```{r eval = FALSE}
src3 <- bb_source(
    name = "Sea Ice Trends and Climatologies from SMMR and SSM/I-SSMIS, Version 3",
    ... details as above...,
    user = "",
    password = "",
    warn_empty_auth = FALSE)
```

When another user wants to use this data source, they simply insert their own credentials, e.g.:

```{r eval = FALSE}
mysrc <- src3
mysrc$user <- "theirusername"
mysrc$password <- "theirpassword"

## or, do the same with the pipe operator and bb_modify_source
mysrc <- src3 %>% bb_modify_source(user = "theirusername", password = "theirpassword")

## then proceed as per usual
cf <- bb_add(cf, mysrc)
```

#### Example 4: an Oceandata source

NASA's [Oceandata](https://oceandata.sci.gsfc.nasa.gov/) system provides access to a range of satellite-derived marine data products. The `bb_oceandata_handler` can be used to download these data. It uses a two-step process: first it makes a query to the Oceancolour data file search tool (https://oceandata.sci.gsfc.nasa.gov/search/file_search.cgi) to find files that match your specified criterion, and then downloads the matching files.

Oceandata uses standardized file naming conventions (see https://oceancolor.gsfc.nasa.gov/docs/format/), so once you know which products you want you can construct a suitable file name pattern to search for. For example, "S*L3m_MO_CHL_chlor_a_9km.nc" would match monthly level-3 mapped chlorophyll data from the SeaWiFS satellite at 9km resolution, in netcdf format. This pattern is passed as the `search` argument to the `bb_handler_oceandata` handler function. Note that the `bb_handler_oceandata` does not need a `source_url` to be specified in the `bb_source` call.

Note that as of 1-Jan-2020, all Oceandata users [require an EarthData login](https://oceancolor.gsfc.nasa.gov/forum/oceancolor/topic_show.pl?tid=11280). You will need to provide your username and password in the `bb_source` call, and you will also need to authorize the 'OB.DAAC Data Access' application in your EarthData account (see 'My Applications' at https://urs.earthdata.nasa.gov/profile).

Here, for the sake of a small example, we'll limit ourselves to a single file ("T20000322000060.L3m_MO_PAR_par_9km.nc", which is photosynthetically-active radiation from the Terra satellite in February 2000):

```{r eval = FALSE}
src4 <- bb_source(
    name = "Oceandata test file",
    id = "oceandata-test",
    description = "Monthly, 9km remote-sensed PAR from the MODIS Terra satellite",
    doc_url = "https://oceancolor.gsfc.nasa.gov/",
    citation = "See https://oceancolor.gsfc.nasa.gov/cms/citations",
    license = "Please cite",
    method = list("bb_handler_oceandata", search="T20000322000060.L3m_MO_PAR_par_9km.nc"),
    user = "your_earthdata_username",
    password = "your_earthdata_password",
    data_group = "Ocean data")

## add this source to a configuration and synchronize it:
cf <- bb_config("c:/temp/data/bbtest") %>% bb_add(src4)
status <- bb_sync(cf)

## and our local copy of this data file:
status$files[[1]]$file
```

## Nuances

Bowerbird hands off the complexities of recursive downloading to the `bb_rget` utility. This allows bowerbird's data source definitions to be fairly lightweight and more robust to changes by the data provider. However, one of the consequences of this approach is that bowerbird actually knows very little about the data files that it maintains, which can be limiting in some respects. It is not generally possible, for example, to provide the user with an indication of overall download progress (progress bar or similar) for a given data source because neither bowerbird nor `bb_rget` actually know in advance how many files are in it or how big they are. Data sources do have a `collection_size` entry, to give the user some indication of the disk space required, but this is only approximate (and must be hand-coded by the data source maintainer). See the 'Reducing download sizes' section below for tips on retrieving only a subset of a large data source.

### Recursive download

- the synchronization process saves files relative to the `local_file_root` directory specified in `bb_config`. The `bb_rget` function saves files into a directory structure that follows the URL structure. For example, calling `bb_rget("http://www.somewhere.org/monkey/banana/dataset.zip")` will save the local file `www.somewhere.org/monkey/banana/dataset.zip`. Thus, `bb_rget` will keep data files from different sources naturally separated into their own directories

Recursion is a powerful tool but will sometimes download much more than you really wanted. There are various methods for restricting the recursion:

- if you want to include/exclude certain files from being downloaded, use the `accept_follow`, `reject_follow`, `accept_download`, and `reject_download` parameters. The `*_follow` parameters control which links will be followed by the recursive spidering process, whereas the `*_download` parameters control which data files will be downloaded

- `no_parent = TRUE` prevents `bb_rget` from ascending to a parent directory during its recursion process, because if it did so it would likely be downloading files that are not part of the data set that we want (this is `TRUE` by default). In some cases, though, you might want the recursion to ascend to a parent directory, and therefore need to override the default setting


#### Other tips and tricks, including resolving recursive download issues

- `no_check_certificate = TRUE` will allow a download from a secure server to proceed even if the server's certificate checks fail. This option might be useful if trying to download files from a server with an expired certificate, but it is clearly a security risk and so should be used with caution

- setting `wait` will cause `bb_rget` to pause for this number of seconds between successive retrievals. This option may help with servers that block multiple successive requests, by introducing a delay between requests

- if `bb_rget` is not behaving as expected, try adding `debug = TRUE`. This gives debugging output from the underlying `libcurl` calls (which is additional to the output obtained by calling `bb_sync(...,verbose=TRUE)`).


### Choosing a data directory

It's up to you where you want your data collection kept, and to provide that location to bowerbird. A common use case for bowerbird is maintaining a central data collection for multiple users, in which case that location is likely to be some sort of networked file share. However, if you are keeping a collection for your own use, you might like to look at https://github.com/r-lib/rappdirs to help find a suitable directory location.


### Post-processing

#### Decompressing files

If the data source delivers compressed files, you will most likely want to decompress them after downloading. The postprocess options `bb_decompress`, `bb_unzip`, etc will do this for you. By default, these *do not* delete the compressed files after decompressing. The reason for this is so that on the next synchronization run, the local (compressed) copy can be compared to the remote compressed copy, and the download can be skipped if nothing has changed. Deleting local compressed files will save space on your file system, but may result in every file being re-downloaded on every synchronization run.

See `help("bb_unzip")` for more information, including usage examples.

#### Deleting unwanted files

The `bb_cleanup` postprocessing option can be used to remove unwanted files after downloading. See See `help("bb_cleanup")`.


### Modifying data sources

#### Authentication

Some data providers require users to log in. The `authentication_note` column in the configuration table should indicate when this is the case, including a reference (e.g. the URL via which an account can be obtained). For these sources, you will need to provide your user name and password, e.g.:

```{r eval = FALSE}
mysrc <- bb_example_sources("CMEMS global gridded SSH reprocessed (1993-ongoing)")
mysrc$user <- "yourusername"
mysrc$password <- "yourpassword"
cf <- bb_add(cf, mysrc)

## or, using the pipe operator
mysrc <- bb_example_sources("CMEMS global gridded SSH reprocessed (1993-ongoing)") %>%
  bb_modify_source(user = "yourusername", password = "yourpassword")
cf <- cf %>% bb_add(mysrc)
```

#### Reducing download sizes

Sometimes you might only want part of a data collection. Perhaps you only want a few years from a long-term collection, or perhaps the data are provided in multiple formats and you only need one. If the data source uses the `bb_handler_rget` method, you can restrict what is downloaded by modifying the arguments passed through the data source's `method` parameter, particularly the `accept_follow`, `reject_follow`, `accept_download`, and `reject_download` options. If you are modifying an existing data source configuration, you most likely want to leave the original method flags intact and just add extra flags.

Say a particular data provider uses predictable file name patterns that include the year information. It would be fairly easy to restrict ourselves to only the 2017 data using the `accept` option. Here we use the `bb_modify_source` helper function to do so:

```{r eval = FALSE}
mysrc <- bb_example_sources("CMEMS global gridded SSH reprocessed (1993-ongoing)") %>%
  bb_modify_source(user = "yourusername", password = "yourpassword",
    method = list(accept_follow = "/2017"))

cf <- cf %>% bb_add(mysrc)
```

Alternatively, for data sources that are arranged in subdirectories, one could replace the `source_url` with one or more that point to the specific subdirectories that are wanted.


### Parallelized sync

If you have many data sources in your configuration, running the sync in parallel is likely to speed the process up considerably (unless your bandwidth is the limiting factor). A logical approach to this would be to split a configuration, with a subset of data sources in each (see `bb_subset`), and run those subsets in parallel. One potential catch to keep in mind would be data sources that hit the same remote data provider. If they overlap overlap in terms of the parts of the remote site that they are mirroring, that might invoke odd behaviour (race conditions, simultaneous downloads of the same file by different parallel processes, etc).


### Data provenance and reproducible research

An aspect of reproducible research is knowing which data were used to perform an analysis, and potentially archiving those data to an appropriate repository. Bowerbird can assist with this: see `vignette("data_provenance")`.


## Developer notes

### Writing new data source handlers

The `bb_handler_rget` R function provides a wrapper around `bb_rget` that should be sufficient for many data sources. However, some data sources need a more elaborate method. Notes will be added here about defining new handler functions, but in the meantime look at `bb_handler_oceandata` and `bb_handler_earthdata`, which provide handlers for [Oceandata](https://oceandata.sci.gsfc.nasa.gov/) and [Earthdata](https://earthdata.nasa.gov/) data sources.

## Data source summary

The `bb_summary` function will produce a HTML or Rmarkdown summary of the data sources contained in a configuration object. If you are maintaining a data collection on behalf of other users, or even just for yourself, it may be useful to keep an up-to-date HTML summary of your repository in an accessible location. Users can refer to this summary to see which data are in the repository and some details about them.

Here is a `bb_summary` of the example data source definitions that are provided as part of the bowerbird package:

```{r echo = FALSE, message = FALSE, warning = FALSE, results = "asis"}
library(bowerbird)
cf <- bb_config("/your/data/root/") %>% bb_add(bb_example_sources())
sf <- bb_summary(cf, format = "rmd", inc_license = FALSE, inc_access_function = FALSE, inc_path = FALSE)
stxt <- readLines(sf)
stxt <- stxt[(grep("Last updated:",stxt)+1):length(stxt)]
stxt <- gsub("^#", "##", stxt) ## push each header level down one
stxt <- gsub("^\\-", "\n-", stxt)
for (k in stxt) cat(k, "\n")
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wget.R
\name{bb_find_wget}
\alias{bb_find_wget}
\title{Find the wget executable}
\usage{
bb_find_wget(install = FALSE, error = TRUE)
}
\arguments{
\item{install}{logical: attempt to install the executable if it is not found? (Windows only)}

\item{error}{logical: if wget is not found, raise an error. If \code{FALSE}, do not raise an error but return NULL}
}
\value{
the path to the wget executable, or (if error is \code{FALSE}) NULL if it was not found
}
\description{
This function will return the path to the wget executable if it can be found on the local system, and optionally install it if it is not found. Installation (if required) currently only works on Windows platforms. The wget.exe executable will be downloaded from https://eternallybored.org/misc/wget/ installed into your appdata directory (typically something like C:/Users/username/AppData/Roaming/)
}
\examples{
\dontrun{
  wget_path <- bb_find_wget()
  wget_path <- bb_find_wget(install=TRUE) ## install (on windows) if needed
}

}
\references{
https://eternallybored.org/misc/wget/
}
\seealso{
\code{\link{bb_install_wget}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zenodo.R
\name{bb_zenodo_source}
\alias{bb_zenodo_source}
\title{Generate a bowerbird data source object for a Zenodo data set}
\usage{
bb_zenodo_source(id)
}
\arguments{
\item{id}{: the ID of the data set}
}
\value{
A tibble containing the data source definition, as would be returned by \code{\link{bb_source}}
}
\description{
Generate a bowerbird data source object for a Zenodo data set
}
\examples{
\dontrun{
  ## generate the source object for the dataset
  ##   'Ichtyological data of Station de biologie des Laurentides 2019'
  src <- bb_zenodo_source(3533328)

  ## download it to a temporary directory
  data_dir <- tempfile()
  dir.create(data_dir)
  res <- bb_get(src, local_file_root = data_dir, verbose = TRUE)
  res$files
}

}
\seealso{
\code{\link{bb_source}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/example_sources.R
\name{bb_example_sources}
\alias{bb_example_sources}
\title{Example bowerbird data sources}
\usage{
bb_example_sources(sources)
}
\arguments{
\item{sources}{character: names or identifiers of one or more sources to return. See Details for the list of example sources and a brief explanation of each}
}
\value{
a tibble with columns as specified by \code{\link{bb_source}}
}
\description{
These example sources are useful as data sources in their own right, but are primarily provided as demonstrations of how to define data sources. See also \code{vignette("bowerbird")} for further examples and discussion.
}
\details{
Example data sources:
\itemize{
  \item "NOAA OI SST V2" - a straightforward data source that requires a simple one-level recursive download
  \item "Australian Election 2016 House of Representatives data" - an example of a recursive download that uses additional criteria to restrict what is downloaded
  \item "CMEMS global gridded SSH reprocessed (1993-ongoing)" - a data source that requires a username and password
  \item "Oceandata SeaWiFS Level-3 mapped monthly 9km chl-a" - an example data source that uses the \code{bb_handler_oceandata} method
  \item "Sea Ice Trends and Climatologies from SMMR and SSM/I-SSMIS, Version 3" - an example data source that uses the \code{bb_handler_earthdata} method
  \item "Bathymetry of Lake Superior" - another example that passes extra flags to the \code{bb_handler_rget} call in order to restrict what is downloaded
}
}
\examples{
## define a configuration and add the 2016 election data source to it
cf <- bb_config("/my/file/root") \%>\% bb_add(
   bb_example_sources("Australian Election 2016 House of Representatives data"))

\dontrun{
  ## synchronize (download) the data
  bb_sync(cf)
}
}
\references{
See the \code{doc_url} and \code{citation} field in each row of the returned tibble for references associated with these particular data sources
}
\seealso{
\code{\link{bb_config}}, \code{\link{bb_handler_rget}}, \code{\link{bb_handler_oceandata}}, \code{\link{bb_handler_earthdata}}, \code{\link{bb_source_us_buildings}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wget.R
\name{bb_handler_wget}
\alias{bb_handler_wget}
\title{Mirror an external data source using the wget utility}
\usage{
bb_handler_wget(...)
}
\arguments{
\item{...}{: parameters passed to \code{\link{bb_wget}}}
}
\value{
TRUE on success
}
\description{
This is a general handler function that is suitable for a range of data sets. This function is not intended to be called directly, but rather is specified as a \code{method} option in \code{\link{bb_source}}.
}
\details{
This handler function makes calls to the \code{wget} utility via the \code{\link{bb_wget}} function. Arguments provided to \code{bb_handler_wget} are passed through to \code{\link{bb_wget}}.
}
\examples{

my_source <- bb_source(
   id="gshhg_coastline",
   name="GSHHG coastline data",
   description="A Global Self-consistent, Hierarchical, High-resolution Geography Database",
   doc_url= "http://www.soest.hawaii.edu/pwessel/gshhg",
   citation="Wessel, P., and W. H. F. Smith, A Global Self-consistent, Hierarchical,
     High-resolution Shoreline Database, J. Geophys. Res., 101, 8741-8743, 1996",
   source_url="ftp://ftp.soest.hawaii.edu/gshhg/*",
   license="LGPL",
   method=list("bb_handler_wget",recursive=TRUE,level=1,accept="*bin*.zip,README.TXT"),
   postprocess=list("bb_unzip"),
   collection_size=0.6)

}
\seealso{
\code{\link{bb_wget}}, \code{\link{bb_source}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/config.R
\name{bb_summary}
\alias{bb_summary}
\title{Produce a summary of a bowerbird configuration}
\usage{
bb_summary(
  config,
  file = tempfile(fileext = ".html"),
  format = "html",
  inc_license = TRUE,
  inc_auth = TRUE,
  inc_size = TRUE,
  inc_access_function = TRUE,
  inc_path = TRUE
)
}
\arguments{
\item{config}{bb_config: a bowerbird configuration (as returned by \code{bb_config})}

\item{file}{string: path to file to write summary to. A temporary file is used by default}

\item{format}{string: produce HTML ("html") or Rmarkdown ("Rmd") file?}

\item{inc_license}{logical: include each source's license and citation details?}

\item{inc_auth}{logical: include information about authentication for each data source (if applicable)?}

\item{inc_size}{logical: include each source's size (disk space) information?}

\item{inc_access_function}{logical: include each source's access function?}

\item{inc_path}{logical: include each source's local file path?}
}
\value{
path to the summary file in HTML or Rmarkdown format
}
\description{
This function produces a summary of a bowerbird configuation in HTML or Rmarkdown format. If you are maintaining a data collection on behalf of other users, or even just for yourself, it may be useful to keep an up-to-date HTML summary of your repository in an accessible location. Users can refer to this summary to see which data are in the repository and some details about them.
}
\examples{
\dontrun{
  cf <- bb_config("/my/file/root") \%>\%
    bb_add(bb_example_sources())
  browseURL(bb_summary(cf))
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/postprocess.R
\name{bb_decompress}
\alias{bb_decompress}
\alias{bb_unzip}
\alias{bb_gunzip}
\alias{bb_bunzip2}
\alias{bb_uncompress}
\alias{bb_untar}
\title{Postprocessing: decompress zip, gz, bz2, tar, Z files and optionally delete the compressed copy}
\usage{
bb_decompress(method, delete = FALSE, ...)

bb_unzip(...)

bb_gunzip(...)

bb_bunzip2(...)

bb_uncompress(...)

bb_untar(...)
}
\arguments{
\item{method}{string: one of "unzip", "gunzip", "bunzip2", "decompress", "untar"}

\item{delete}{logical: delete the zip files after extracting their contents?}

\item{...}{: extra parameters passed automatically by \code{bb_sync}}
}
\value{
list with components status (\code{TRUE} on success), \code{files} (character vector of paths to extracted files), and \code{deleted_files} (character vector of paths of files that were deleted)
}
\description{
Functions for decompressing files after downloading. These functions are not intended to be called directly, but rather are specified as a \code{postprocess} option in \code{\link{bb_source}}. \code{bb_unzip}, \code{bb_untar}, \code{bb_gunzip}, \code{bb_bunzip2}, and \code{bb_uncompress} are convenience wrappers around \code{bb_decompress} that specify the method.
}
\details{
Tar files can be compressed (i.e. file extensions .tar, .tgz, .tar.gz, .tar.bz2, or .tar.xz). Support for tar files may depend on your platform (see \code{\link{untar}}).

If the data source delivers compressed files, you will most likely want to decompress them after downloading. These functions will do this for you. By default, these do not delete the compressed files after decompressing. The reason for this is so that on the next synchronization run, the local (compressed) copy can be compared to the remote compressed copy, and the download can be skipped if nothing has changed. Deleting local compressed files will save space on your file system, but may result in every file being re-downloaded on every synchronization run.
}
\examples{
\dontrun{
  ## decompress .zip files after synchronization but keep zip files intact
  my_source <- bb_source(..., postprocess = list("bb_unzip"))

  ## decompress .zip files after synchronization and delete zip files
  my_source <- bb_source(..., postprocess = list(list("bb_unzip", delete = TRUE)))
}

}
\seealso{
\code{\link{bb_source}}, \code{\link{bb_config}}, \code{\link{bb_cleanup}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rget.R
\name{bb_rget}
\alias{bb_rget}
\alias{bb_rget_default_downloads}
\title{A recursive download utility}
\usage{
bb_rget(
  url,
  level = 0,
  wait = 0,
  accept_follow = c("(/|\\\\.html?)$"),
  reject_follow = character(),
  accept_download = bb_rget_default_downloads(),
  accept_download_extra = character(),
  reject_download = character(),
  user,
  password,
  clobber = 1,
  no_parent = TRUE,
  no_check_certificate = FALSE,
  relative = FALSE,
  remote_time = TRUE,
  verbose = FALSE,
  show_progress = verbose,
  debug = FALSE,
  dry_run = FALSE,
  stop_on_download_error = FALSE,
  force_local_filename,
  use_url_directory = TRUE,
  no_host = FALSE,
  cut_dirs = 0L,
  link_css = "a",
  curl_opts
)

bb_rget_default_downloads()
}
\arguments{
\item{url}{string: the URL to retrieve}

\item{level}{integer >=0: recursively download to this maximum depth level. Specify 0 for no recursion}

\item{wait}{numeric >=0: wait this number of seconds between successive retrievals. This option may help with servers that block users making too many requests in a short period of time}

\item{accept_follow}{character: character vector with one or more entries. Each entry specifies a regular expression that is applied to the complete URL. URLs matching all entries will be followed during the spidering process. Note that the first URL (provided via the \code{url} parameter) will always be visited, unless it matches the download criteria}

\item{reject_follow}{character: as for \code{accept_follow}, but specifying URL regular expressions to reject}

\item{accept_download}{character: character vector with one or more entries. Each entry specifies a regular expression that is applied to the complete URL. URLs that match all entries will be accepted for download. By default the \code{accept_download} parameter is that returned by \code{bb_rget_default_downloads}: use \code{bb_rget_default_downloads()} to see what that is}

\item{accept_download_extra}{character: character vector with one or more entries. If provided, URLs will be accepted for download if they match all entries in \code{accept_download} OR all entries in \code{accept_download_extra}. This is a convenient method to add one or more extra download types, without needing to re-specify the defaults in \code{accept_download}}

\item{reject_download}{character: as for \code{accept_regex}, but specifying URL regular expressions to reject}

\item{user}{string: username used to authenticate to the remote server}

\item{password}{string: password used to authenticate to the remote server}

\item{clobber}{numeric: 0=do not overwrite existing files, 1=overwrite if the remote file is newer than the local copy, 2=always overwrite existing files}

\item{no_parent}{logical: if \code{TRUE}, do not ever ascend to the parent directory when retrieving recursively. This is \code{TRUE} by default, bacause it guarantees that only the files below a certain hierarchy will be downloaded. Note that this check only applies to links on the same host as the starting \code{url}. If that URL links to files on another host, those links will be followed (unless \code{relative = TRUE})}

\item{no_check_certificate}{logical: if \code{TRUE}, don't check the server certificate against the available certificate authorities. Also don't require the URL host name to match the common name presented by the certificate. This option might be useful if trying to download files from a server with an expired certificate, but it is clearly a security risk and so should be used with caution}

\item{relative}{logical: if \code{TRUE}, only follow relative links. This can be useful for restricting what is downloaded in recursive mode}

\item{remote_time}{logical: if \code{TRUE}, attempt to set the local file's time to that of the remote file}

\item{verbose}{logical: print trace output?}

\item{show_progress}{logical: if \code{TRUE}, show download progress}

\item{debug}{logical: if \code{TRUE}, will print additional debugging information. If bb_rget is not behaving as expected, try setting this to \code{TRUE}}

\item{dry_run}{logical: if \code{TRUE}, spider the remote site and work out which files would be downloaded, but don't download them}

\item{stop_on_download_error}{logical: if \code{TRUE}, the download process will stop if any file download fails. If \code{FALSE}, the process will issue a warning and continue to the next file to download}

\item{force_local_filename}{string: if provided, then the \code{url} will be treated as a single URL (no recursion will be conducted). It will be downloaded to a file with this name, in a local directory determined by the \code{url}}

\item{use_url_directory}{logical: if \code{TRUE}, files will be saved into a local directory that follows the URL structure (e.g. files from \code{http://some.where/place} will be saved into directory \code{some.where/place}). If \code{FALSE}, files will be saved into the current directory}

\item{no_host}{logical: if \code{use_url_directory = TRUE}, specifying \code{no_host = TRUE} will remove the host name from the directory (e.g. files from files from \code{http://some.where/place} will be saved into directory \code{place})}

\item{cut_dirs}{integer: if \code{use_url_directory = TRUE}, specifying \code{cut_dirs} will remove this many directory levels from the path of the local directory where files will be saved (e.g. if \code{cut_dirs = 2}, files from \code{http://some.where/place/baa/haa} will be saved into directory \code{some.where/haa}. if \code{cut_dirs = 1} and \code{no_host = TRUE}, files from \code{http://some.where/place/baa/haa} will be saved into directory \code{baa/haa})}

\item{link_css}{string: css selector that identifies links (passed as the \code{css} parameter to \code{\link[rvest]{html_nodes}}). Note that link elements must have an \code{href} attribute}

\item{curl_opts}{named list: options to use with \code{curl} downloads, passed to the \code{.list} parameter of \code{curl::new_handle}}
}
\value{
a list with components 'ok' (TRUE/FALSE), 'files', and 'message' (error or other messages)
}
\description{
This function provides similar, but simplified, functionality to the the command-line \code{wget} utility. It is based on the \code{rvest} package.
}
\details{
NOTE: this is still somewhat experimental.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rget.R
\name{bb_handler_rget}
\alias{bb_handler_rget}
\title{Mirror an external data source using bowerbird's bb_rget utility}
\usage{
bb_handler_rget(...)
}
\arguments{
\item{...}{: parameters passed to \code{\link{bb_rget}}}
}
\value{
TRUE on success
}
\description{
This is a general handler function that is suitable for a range of data sets. This function is not intended to be called directly, but rather is specified as a \code{method} option in \code{\link{bb_source}}.
}
\details{
This handler function makes calls to the \code{\link{bb_rget}} function. Arguments provided to \code{bb_handler_rget} are passed through to \code{\link{bb_rget}}.
}
\examples{
my_source <- bb_source(
   name = "Australian Election 2016 House of Representatives data",
   id = "aus-election-house-2016",
   description = "House of Representatives results from the 2016 Australian election.",
   doc_url = "http://results.aec.gov.au/",
   citation = "Copyright Commonwealth of Australia 2017. As far as practicable, material for
               which the copyright is owned by a third party will be clearly labelled. The
               AEC has made all reasonable efforts to ensure that this material has been
               reproduced on this website with the full consent of the copyright owners.",
   source_url = "http://results.aec.gov.au/20499/Website/HouseDownloadsMenu-20499-Csv.htm",
   license = "CC-BY",
   method = list("bb_handler_rget", level = 1, accept_download = "csv$"),
   collection_size = 0.01)

my_data_dir <- tempdir()
cf <- bb_config(my_data_dir)
cf <- bb_add(cf, my_source)

\dontrun{
bb_sync(cf, verbose = TRUE)
}

}
\seealso{
\code{\link{bb_rget}}, \code{\link{bb_source}}, \code{\link{bb_sync}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/config.R
\name{bb_add}
\alias{bb_add}
\title{Add new data sources to a bowerbird configuration}
\usage{
bb_add(config, source)
}
\arguments{
\item{config}{bb_config: a bowerbird configuration (as returned by \code{bb_config})}

\item{source}{data.frame: one or more data source definitions, as returned by \code{bb_source}, to add to the configuration}
}
\value{
configuration object
}
\description{
Add new data sources to a bowerbird configuration
}
\examples{
\dontrun{
  cf <- bb_config("/my/file/root") \%>\%
    bb_add(bb_example_sources())
}
}
\seealso{
\code{\link{bb_source}}, \code{\link{bb_config}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wget.R
\name{bb_install_wget}
\alias{bb_install_wget}
\title{Install wget}
\usage{
bb_install_wget(force = FALSE, use_appdata_dir = FALSE)
}
\arguments{
\item{force}{logical: force reinstallation if wget already exists}

\item{use_appdata_dir}{logical: by default, \code{bb_install_wget} will install wget into a temporary directory, which does not persist between R sessions. If you want a persistent installation, specify \code{use_appdata_dir=TRUE} to install wget into your appdata directory (on Windows, typically something like C:/Users/username/AppData/Roaming/)}
}
\value{
the path to the installed executable
}
\description{
This is a helper function to install wget. Currently it only works on Windows platforms. The wget.exe executable will be downloaded from https://eternallybored.org/misc/wget/ and saved to either a temporary directory or your user appdata directory (see the \code{use_appdata_dir} parameter).
}
\examples{
\dontrun{
  bb_install_wget()

  ## confirm that it worked:
  bb_wget("help")
}

}
\references{
https://eternallybored.org/misc/wget/
}
\seealso{
\code{\link{bb_find_wget}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/source.R
\name{bb_modify_source}
\alias{bb_modify_source}
\title{Modify a data source}
\usage{
bb_modify_source(src, ...)
}
\arguments{
\item{src}{data.frame or tibble: a single-row data source (as returned by \code{bb_source})}

\item{...}{: parameters as for \code{bb_source}}
}
\value{
as for \code{bb_source}: a tibble with columns as per the \code{bb_source} function arguments (excluding \code{warn_empty_auth})
}
\description{
This is a helper function designed to make it easier to modify an already-defined data source. Generally, parameters passed here will replace existing entries in \code{src} if they exist, or will be added if not. The \code{method} and \code{postprocess} parameters are slightly different: see Details, below.
}
\details{
With the exception of the \code{method} and \code{postprocess} parameters, any parameter provided here will entirely replace its equivalent in the \code{src} object. Pass a new value of \code{NULL} to remove an existing parameter.

The \code{method} and \code{postprocess} parameters are lists, and modification for these takes place at the list-element level: any element of the new list will replace its equivalent element in the list in src. If the src list does not contain that element, it will be added. To illustrate, say that we have created a data source with:

\code{src <- bb_source(method=list("bb_handler_rget", parm1 = value1, parm2 = value2), ...)}

Calling

\code{bb_modify_source(src, method = list(parm1 = newvalue1))}

will result in a new \code{method} value of \code{list("bb_handler_rget", parm1 = newvalue1, parm2 = value2)}

Modifying \code{postprocess} elements is similar. Note that it is not currently possible to entirely remove a postprocess component using this function. If you need to do so, you'll need to do it manually.
}
\examples{

## this pre-defined source requires a username and password
src <- bb_example_sources(
          "Sea Ice Trends and Climatologies from SMMR and SSM/I-SSMIS, Version 3")

## add username and password
src <- bb_modify_source(src,user="myusername",password="mypassword")

## or using the pipe operator
src <- bb_example_sources(
          "Sea Ice Trends and Climatologies from SMMR and SSM/I-SSMIS, Version 3") \%>\%
  bb_modify_source(user="myusername",password="mypassword")

## remove the existing "data_group" component
src \%>\% bb_modify_source(data_group=NULL)

## change just the 'level' setting of an existing method definition
src \%>\% bb_modify_source(method=list(level=3))

## remove the 'level' component of an existing method definition
src \%>\% bb_modify_source(method=list(level=NULL))

}
\seealso{
\code{\link{bb_source}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aadc.R
\name{bb_aadc_source}
\alias{bb_aadc_source}
\title{Generate a bowerbird data source object for an Australian Antarctic Data Centre data set}
\usage{
bb_aadc_source(metadata_id)
}
\arguments{
\item{metadata_id}{string: the metadata ID of the data set. Browse the AADC's collection at https://data.aad.gov.au/metadata/records/ to find the relevant \code{metadata_id}}
}
\value{
A tibble containing the data source definition, as would be returned by \code{\link{bb_source}}
}
\description{
Generate a bowerbird data source object for an Australian Antarctic Data Centre data set
}
\examples{
\dontrun{
  ## generate the source def for the "AADC-00009" dataset
  ##  (Antarctic Fur Seal Populations on Heard Island, Summer 1987-1988)
  src <- bb_aadc_source("AADC-00009")

  ## download it to a temporary directory
  data_dir <- tempfile()
  dir.create(data_dir)
  res <- bb_get(src, local_file_root = data_dir, verbose = TRUE)
  res$files
}

}
\seealso{
\code{\link{bb_source}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wget.R
\name{bb_wget}
\alias{bb_wget}
\title{Make a wget call}
\usage{
bb_wget(
  url,
  recursive = TRUE,
  level = 1,
  wait = 0,
  accept,
  reject,
  accept_regex,
  reject_regex,
  exclude_directories,
  restrict_file_names,
  progress,
  user,
  password,
  output_file,
  robots_off = FALSE,
  timestamping = FALSE,
  no_if_modified_since = FALSE,
  no_clobber = FALSE,
  no_parent = TRUE,
  no_check_certificate = FALSE,
  relative = FALSE,
  adjust_extension = FALSE,
  retr_symlinks = FALSE,
  extra_flags = character(),
  verbose = FALSE,
  capture_stdout = FALSE,
  quiet = FALSE,
  debug = FALSE
)
}
\arguments{
\item{url}{string: the URL to retrieve}

\item{recursive}{logical: if true, turn on recursive retrieving}

\item{level}{integer >=0: recursively download to this maximum depth level. Only applicable if \code{recursive=TRUE}. Specify 0 for infinite recursion. See \url{https://www.gnu.org/software/wget/manual/wget.html#Recursive-Download} for more information about wget's recursive downloading}

\item{wait}{numeric >=0: wait this number of seconds between successive retrievals. This option may help with servers that block multiple successive requests, by introducing a delay between requests}

\item{accept}{character: character vector with one or more entries. Each entry specifies a comma-separated list of filename suffixes or patterns to accept. Note that if any of the wildcard characters '*', '?', '[', or ']' appear in an element of accept, it will be treated as a filename pattern, rather than a filename suffix. In this case, you have to enclose the pattern in quotes, for example \code{accept="\"*.csv\""}}

\item{reject}{character: as for \code{accept}, but specifying filename suffixes or patterns to reject}

\item{accept_regex}{character: character vector with one or more entries. Each entry provides a regular expression that is applied to the complete URL. Matching URLs will be accepted for download}

\item{reject_regex}{character: as for \code{accept_regex}, but specifying regular expressions to reject}

\item{exclude_directories}{character: character vector with one or more entries. Each entry specifies a comma-separated list of directories you wish to exclude from download. Elements may contain wildcards}

\item{restrict_file_names}{character: vector of one of more strings from the set "unix", "windows", "nocontrol", "ascii", "lowercase", and "uppercase". See \url{https://www.gnu.org/software/wget/manual/wget.html#index-Windows-file-names} for more information on this parameter. \code{bb_config} sets this to "windows" by default: if you are downloading files from a server with a port (http://somewhere.org:1234/) Unix will allow the ":" as part of directory/file names, but Windows will not (the ":" will be replaced by "+"). Specifying \code{restrict_file_names="windows"} causes Windows-style file naming to be used}

\item{progress}{string: the type of progress indicator you wish to use. Legal indicators are "dot" and "bar". "dot" prints progress with dots, with each dot representing a fixed amount of downloaded data. The style can be adjusted: "dot:mega" will show 64K per dot and 3M per line; "dot:giga" shows 1M per dot and 32M per line. See \url{https://www.gnu.org/software/wget/manual/wget.html#index-dot-style} for more information}

\item{user}{string: username used to authenticate to the remote server}

\item{password}{string: password used to authenticate to the remote server}

\item{output_file}{string: save wget's output messages to this file}

\item{robots_off}{logical: by default wget considers itself to be a robot, and therefore won't recurse into areas of a site that are excluded to robots. This can cause problems with servers that exclude robots (accidentally or deliberately) from parts of their sites containing data that we want to retrieve. Setting \code{robots_off=TRUE} will add a "-e robots=off" flag, which instructs wget to behave as a human user, not a robot. See \url{https://www.gnu.org/software/wget/manual/wget.html#Robot-Exclusion} for more information about robot exclusion}

\item{timestamping}{logical: if \code{TRUE}, don't re-retrieve a remote file unless it is newer than the local copy (or there is no local copy)}

\item{no_if_modified_since}{logical: applies when retrieving recursively with timestamping (i.e. only downloading files that have changed since last download, which is achieved using \code{bb_config(...,clobber=1)}). The default method for timestamping is to issue an "If-Modified-Since" header on the request, which instructs the remote server not to return the file if it has not changed since the specified date. Some servers do not support this header. In these cases, trying using \code{no_if_modified_since=TRUE}, which will instead send a preliminary HEAD request to ascertain the date of the remote file}

\item{no_clobber}{logical: if \code{TRUE}, skip downloads that would overwrite existing local files}

\item{no_parent}{logical: if \code{TRUE}, do not ever ascend to the parent directory when retrieving recursively. This is \code{TRUE} by default, bacause it guarantees that only the files below a certain hierarchy will be downloaded}

\item{no_check_certificate}{logical: if \code{TRUE}, don't check the server certificate against the available certificate authorities. Also don't require the URL host name to match the common name presented by the certificate. This option might be useful if trying to download files from a server with an expired certificate, but it is clearly a security risk and so should be used with caution}

\item{relative}{logical: if \code{TRUE}, only follow relative links. This can sometimes be useful for restricting what is downloaded in recursive mode}

\item{adjust_extension}{logical: if a file of type 'application/xhtml+xml' or 'text/html' is downloaded and the URL does not end with .htm or .html, this option will cause the suffix '.html' to be appended to the local filename. This can be useful when mirroring a remote site that has file URLs that conflict with directories (e.g. http://somewhere.org/this/page which has further content below it, say at http://somewhere.org/this/page/more. If "somewhere.org/this/page" is saved as a file with that name, that name can't also be used as the local directory name in which to store the lower-level content. Setting \code{adjust_extension=TRUE} will cause the page to be saved as "somewhere.org/this/page.html", thus resolving the conflict}

\item{retr_symlinks}{logical: if \code{TRUE}, follow symbolic links during recursive download. Note that this will only follow symlinks to files, NOT to directories}

\item{extra_flags}{character: character vector of additional command-line flags to pass to wget}

\item{verbose}{logical: print trace output?}

\item{capture_stdout}{logical: if \code{TRUE}, return 'stdout' and 'stderr' output in the returned object (see exec_internal from the sys package). Otherwise send these outputs to the console}

\item{quiet}{logical: if \code{TRUE}, suppress wget's output}

\item{debug}{logical: if \code{TRUE}, wget will print lots of debugging information. If wget is not behaving as expected, try setting this to \code{TRUE}}
}
\value{
the result of the system call (or if \code{bb_wget("--help")} was called, a message will be issued). The returned object will have components 'status' and (if \code{capture_stdout} was \code{TRUE}) 'stdout' and 'stderr'
}
\description{
This function is an R wrapper to the command-line \code{wget} utility, which is called using either the \code{exec_wait} or the \code{exec_internal} function from the sys package. Almost all of the parameters to \code{bb_wget} are translated into command-line flags to \code{wget}. Call \code{bb_wget("help")} to get more information about wget's command line flags. If required, command-line flags without equivalent \code{bb_wget} function parameters can be passed via the \code{extra_flags} parameter.
}
\examples{
\dontrun{
  ## get help about wget command line parameters
  bb_wget("help")
}

}
\seealso{
\code{\link{bb_install_wget}}, \code{\link{bb_find_wget}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/source.R
\name{bb_source}
\alias{bb_source}
\title{Define a data source}
\usage{
bb_source(
  id,
  name,
  description = NA_character_,
  doc_url,
  source_url,
  citation,
  license,
  comment = NA_character_,
  method,
  postprocess,
  authentication_note = NA_character_,
  user = NA_character_,
  password = NA_character_,
  access_function = NA_character_,
  data_group = NA_character_,
  collection_size = NA,
  warn_empty_auth = TRUE
)
}
\arguments{
\item{id}{string: (required) a unique identifier of the data source. If the data source has a DOI, use that. Otherwise, if the original data provider has an identifier for this dataset, that is probably a good choice here (include the data version number if there is one). The ID should be something that changes when the data set changes (is updated). A DOI is ideal for this}

\item{name}{string: (required) a unique name for the data source. This should be a human-readable but still concise name}

\item{description}{string: a plain-language description of the data source, provided so that users can get an idea of what the data source contains (for full details they can consult the \code{doc_url} link)}

\item{doc_url}{string: (required) URL to the metadata record or other documentation of the data source}

\item{source_url}{character vector: one or more source URLs. Required for \code{bb_handler_rget}, although some \code{method} functions might not require one}

\item{citation}{string: (required) details of the citation for the data source}

\item{license}{string: (required) description of the license. For standard licenses (e.g. creative commons) include the license descriptor ("CC-BY", etc)}

\item{comment}{string: comments about the data source. If only part of the original data collection is mirrored, mention that here}

\item{method}{list (required): a list object that defines the function used to synchronize this data source. The first element of the list is the function name (as a string or function). Additional list elements can be used to specify additional parameters to pass to that function. Note that \code{bb_sync} automatically passes the data repository configuration object as the first parameter to the method handler function. If the handler function uses bb_rget (e.g. \code{bb_handler_rget}), these extra parameters are passed through to the \code{bb_rget} function}

\item{postprocess}{list: each element of \code{postprocess} defines a postprocessing step to be run after the main synchronization has happened. Each element of this list can be a function or string function name, or a list in the style of \code{list(fun,arg1=val1,arg2=val2)} where \code{fun} is the function to be called and \code{arg1} and \code{arg2} are additional parameters to pass to that function}

\item{authentication_note}{string: if authentication is required in order to access this data source, make a note of the process (include a URL to the registration page, if possible)}

\item{user}{string: username, if required}

\item{password}{string: password, if required}

\item{access_function}{string: can be used to suggest to users an appropriate function to read these data files. Provide the name of an R function or even a code snippet}

\item{data_group}{string: the name of the group to which this data source belongs. Useful for arranging sources in terms of thematic areas}

\item{collection_size}{numeric: approximate disk space (in GB) used by the data collection, if known. If the data are supplied as compressed files, this size should reflect the disk space used after decompression. If the data_source definition contains multiple source_url entries, this size should reflect the overall disk space used by all combined}

\item{warn_empty_auth}{logical: if \code{TRUE}, issue a warning if the data source requires authentication (authentication_note is not NA) but user and password have not been provided. Set this to \code{FALSE} if you are defining a data source for others to use with their own credentials: they will typically call your data source constructor and then modify the \code{user} and \code{password} components}
}
\value{
a tibble with columns as per the function arguments (excluding \code{warn_empty_auth})
}
\description{
This function is used to define a data source, which can then be added to a bowerbird data repository configuration. Passing the configuration object to \code{bb_sync} will trigger a download of all of the data sources in that configuration.
}
\details{
The \code{method} parameter defines the handler function used to synchronize this data source, and any extra parameters that need to be passed to it.

Parameters marked as "required" are the minimal set needed to define a data source. Other parameters are either not relevant to all data sources (e.g. \code{postprocess}, \code{user}, \code{password}) or provide metadata to users that is not strictly necessary to allow the data source to be synchronized (e.g. \code{description}, \code{access_function}, \code{data_group}). Note that three of the "required" parameters (namely \code{citation}, \code{license}, and \code{doc_url}) are not strictly needed by the synchronization code, but are treated as "required" because of their fundamental importance to reproducible science.

See \code{vignette("bowerbird")} for more examples and discussion of defining data sources.
}
\examples{

## a minimal definition for the GSHHG coastline data set:

my_source <- bb_source(
   id = "gshhg_coastline",
   name = "GSHHG coastline data",
   doc_url = "http://www.soest.hawaii.edu/pwessel/gshhg",
   citation = "Wessel, P., and W. H. F. Smith, A Global Self-consistent, Hierarchical,
     High-resolution Shoreline Database, J. Geophys. Res., 101, 8741-8743, 1996",
   source_url = "ftp://ftp.soest.hawaii.edu/gshhg/",
   license = "LGPL",
   method = list("bb_handler_rget",level = 1, accept_download = "README|bin.*\\\\.zip$"))

## a more complete definition, which unzips the files after downloading and also
##  provides an indication of the size of the dataset

my_source <- bb_source(
   id = "gshhg_coastline",
   name = "GSHHG coastline data",
   description = "A Global Self-consistent, Hierarchical, High-resolution Geography Database",
   doc_url = "http://www.soest.hawaii.edu/pwessel/gshhg",
   citation = "Wessel, P., and W. H. F. Smith, A Global Self-consistent, Hierarchical,
     High-resolution Shoreline Database, J. Geophys. Res., 101, 8741-8743, 1996",
   source_url = "ftp://ftp.soest.hawaii.edu/gshhg/*",
   license = "LGPL",
   method = list("bb_handler_rget", level = 1, accept_download = "README|bin.*\\\\.zip$"),
   postprocess = list("bb_unzip"),
   collection_size = 0.6)

## define a data repository configuration
cf <- bb_config("/my/repo/root")

## add this source to the repository
cf <- bb_add(cf, my_source)

\dontrun{
## sync the repo
bb_sync(cf)
}

}
\seealso{
\code{\link{bb_config}}, \code{\link{bb_sync}}, \code{vignette("bowerbird")}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aws_s3.R
\name{bb_handler_aws_s3}
\alias{bb_handler_aws_s3}
\title{Handler for public AWS S3 data sources}
\usage{
bb_handler_aws_s3(...)
}
\arguments{
\item{...}{: parameters, see Description}
}
\value{
A tibble with columns \code{ok}, \code{files}, \code{message}
}
\description{
This is a handler function to be used with AWS S3 data providers. This function is not intended to be called directly, but rather is specified as a \code{method} option in \code{\link{bb_source}}. Note that this currently only works with public data sources that are accessible without an S3 key.
The method arguments accepted by \code{bb_handler_aws_s3} are currently:
\itemize{
  \item "bucket" string: name of the bucket (defaults to "")
  \item "base_url" string: as for \code{\link[aws.s3]{s3HTTP}}
  \item "region" string: as for \code{\link[aws.s3]{s3HTTP}}
  \item "use_https" logical: as for \code{\link[aws.s3]{s3HTTP}}
  \item "prefix" string: as for \code{\link[aws.s3]{get_bucket}}; only keys in the bucket that begin with the specified prefix will be processed
  \item and other parameters passed to the \code{\link{bb_rget}} function, including "accept_download", "accept_download_extra", "reject_download"
}
Note that the "prefix", "accept_download", "accept_download_extra", "reject_download" parameters can be used to restrict which files are downloaded from the bucket.
}
\examples{
\dontrun{
  ## an example AWS S3 data source
  src <- bb_source(
           name = "SILO climate data",
           id = "silo-open-data",
           description = "Australian climate data from 1889 to yesterday.
                          This source includes a single example monthly rainfall data file.
                          Adjust the 'accept_download' parameter to change this.",
           doc_url = "https://www.longpaddock.qld.gov.au/silo/gridded-data/",
           citation = "SILO datasets are constructed by the Queensland Government using
                       observational data provided by the Australian Bureau of Meteorology
                       and are available under the Creative Commons Attribution 4.0 license",
           license = "CC-BY 4.0",
           method = list("bb_handler_aws_s3", region = "silo-open-data.s3",
                         base_url = "amazonaws.com", prefix = "annual/monthly_rain/",
                         accept_download = "2005\\\\.monthly_rain\\\\.nc$"),
           comment = "The unusual specification of region and base_url is a workaround for
                      an aws.s3 issue, see https://github.com/cloudyr/aws.s3/issues/318",
           postprocess = NULL,
           collection_size = 0.02,
           data_group = "Climate")
   temp_root <- tempdir()
   status <- bb_get(src, local_file_root = temp_root, verbose = TRUE)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/earthdata_handler.R
\name{bb_handler_earthdata}
\alias{bb_handler_earthdata}
\title{Handler for data sets from Earthdata providers}
\usage{
bb_handler_earthdata(...)
}
\arguments{
\item{...}{: parameters passed to \code{\link{bb_rget}}}
}
\value{
TRUE on success
}
\description{
This is a handler function to be used with data sets from NASA's Earthdata system. This function is not intended to be called directly, but rather is specified as a \code{method} option in \code{\link{bb_source}}.
}
\details{
This function uses \code{\link{bb_rget}}, and so data sources using this function will need to provide appropriate \code{\link{bb_rget}} parameters.
}
\examples{
\dontrun{

## note that the full version of this data source is provided as part of bb_example_data_sources()

my_source <- bb_source(
  name = "Sea Ice Trends and Climatologies from SMMR and SSM/I-SSMIS, Version 3",
  id = "10.5067/EYICLBOAAJOU",
  description = "NSIDC provides this data set ... [truncated; see bb_example_data_sources()]",
  doc_url = "https://nsidc.org/data/NSIDC-0192/versions/3",
  citation = "Stroeve J, Meier WN (2018) ... [truncated; see bb_example_data_sources()]",
  source_url = "https://daacdata.apps.nsidc.org/pub/DATASETS/nsidc0192_seaice_trends_climo_v3/",
  license = "Please cite, see http://nsidc.org/about/use_copyright.html",
  authentication_note = "Requires Earthdata login, see https://urs.earthdata.nasa.gov/.
    Note that you will also need to authorize the application 'nsidc-daacdata'
    (see 'My Applications' at https://urs.earthdata.nasa.gov/profile)",
  method = list("bb_handler_earthdata", recursive = TRUE, level = 4, no_parent = TRUE,
                relative = TRUE),
  user = "your_earthdata_username",
  password = "your_earthdata_password",
  collection_size = 0.02)
}

}
\references{
https://wiki.earthdata.nasa.gov/display/EL/How+To+Register+With+Earthdata+Login
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/config.R
\name{bb_subset}
\alias{bb_subset}
\title{Keep only selected data_sources in a bowerbird configuration}
\usage{
bb_subset(config, idx)
}
\arguments{
\item{config}{bb_config: a bowerbird configuration (as returned by \code{bb_config})}

\item{idx}{logical or numeric: index vector of data_source rows to retain}
}
\value{
configuration object
}
\description{
Keep only selected data_sources in a bowerbird configuration
}
\examples{
\dontrun{
  cf <- bb_config("/my/file/root") \%>\%
    bb_add(bb_example_sources()) \%>\%
    bb_subset(1:2)
}
}
\seealso{
\code{\link{bb_source}}, \code{\link{bb_config}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/provenance.R
\name{bb_fingerprint}
\alias{bb_fingerprint}
\title{Fingerprint the files associated with a data source}
\usage{
bb_fingerprint(config, hash = "sha1")
}
\arguments{
\item{config}{bb_config: configuration as returned by \code{\link{bb_config}}}

\item{hash}{string: algorithm to use to calculate file hashes: "md5", "sha1", or "none". Note that file hashing can be slow for large file collections}
}
\value{
a tibble with columns:
\itemize{
  \item filename - the full path and filename of the file
  \item data_source_id - the identifier of the associated data source (as per the \code{id} argument to \code{bb_source})
  \item size - the file size
  \item last_modified - last modified date of the file
  \item hash - the hash of the file (unless \code{hash="none"} was specified)
}
}
\description{
The \code{bb_fingerprint} function, given a data repository configuration, will return the timestamp of download and hashes of all files associated with its data sources. This is intended as a general helper for tracking data provenance: for all of these files, we have information on where they came from (the data source ID), when they were downloaded, and a hash so that later versions of those files can be compared to detect changes. See also \code{vignette("data_provenance")}.
}
\examples{
\dontrun{
  cf <- bb_config("/my/file/root") \%>\%
    bb_add(bb_example_sources())
  bb_fingerprint(cf)
}

}
\seealso{
\code{vignette("data_provenance")}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.R
\name{bb_get}
\alias{bb_get}
\title{Convenience function to define and synchronize a bowerbird data collection}
\usage{
bb_get(
  data_sources,
  local_file_root,
  clobber = 1,
  http_proxy = NULL,
  ftp_proxy = NULL,
  create_root = FALSE,
  verbose = FALSE,
  confirm_downloads_larger_than = 0.1,
  dry_run = FALSE,
  ...
)
}
\arguments{
\item{data_sources}{tibble: one or more data sources to download, as returned by e.g. \code{bb_example_sources}}

\item{local_file_root}{string: location of data repository on local file system}

\item{clobber}{numeric: 0=do not overwrite existing files, 1=overwrite if the remote file is newer than the local copy, 2=always overwrite existing files}

\item{http_proxy}{string: URL of HTTP proxy to use e.g. 'http://your.proxy:8080' (NULL for no proxy)}

\item{ftp_proxy}{string: URL of FTP proxy to use e.g. 'http://your.proxy:21' (NULL for no proxy)}

\item{create_root}{logical: should the data root directory be created if it does not exist? If this is \code{FALSE} (default) and the data root directory does not exist, an error will be generated}

\item{verbose}{logical: if \code{TRUE}, provide additional progress output}

\item{confirm_downloads_larger_than}{numeric or NULL: if non-negative, \code{bb_sync} will ask the user for confirmation to download any data source of size greater than this number (in GB). A value of zero will trigger confirmation on every data source. A negative or NULL value will not prompt for confirmation. Note that this only applies when R is being used interactively. The expected download size is taken from the \code{collection_size} parameter of the data source, and so its accuracy is dependent on the accuracy of the data source definition}

\item{dry_run}{logical: if \code{TRUE}, \code{bb_sync} will do a dry run of the synchronization process without actually downloading files}

\item{...}{: additional parameters passed through to \code{bb_config} or \code{bb_sync}}
}
\value{
a tibble, as for \code{\link{bb_sync}}
}
\description{
This is a convenience function that provides a shorthand method for synchronizing a small number of data sources. The call \code{bb_get(...)} is roughly equivalent to \code{bb_sync(bb_add(bb_config(...), ...), ...)} (don't take the dots literally here, they are just indicating argument placeholders).
}
\details{
Note that the \code{local_file_root} directory must exist or \code{create_root=TRUE} must be passed.
}
\examples{
\dontrun{
  my_source <- bb_example_sources("Australian Election 2016 House of Representatives data")
  status <- bb_get(local_file_root = tempdir(), data_sources = my_source, verbose = TRUE)

  ## the files that have been downloaded:
  status$files[[1]]

  ## Define a new source: Geelong bicycle paths from data.gov.au
  my_source <- bb_source(
    name = "Bike Paths - Greater Geelong",
    id = "http://data.gov.au/dataset/7af9cf59-a4ea-47b2-8652-5e5eeed19611",
    doc_url = "https://data.gov.au/dataset/geelong-bike-paths",
    citation = "See https://data.gov.au/dataset/geelong-bike-paths",
    source_url = "https://data.gov.au/dataset/7af9cf59-a4ea-47b2-8652-5e5eeed19611",
    license = "CC-BY",
    method = list("bb_handler_rget", accept_download = "\\\\.zip$", level = 1),
    postprocess = list("bb_unzip"))

  ## get the data
  status <- bb_get(data_sources = my_source, local_file_root = tempdir(), verbose = TRUE)

  ## find the .shp file amongst the files, and plot it
  shpfile <- status$files[[1]]$file[grepl("shp$", status$files[[1]]$file)]
  library(sf)
  bx <- read_st(shpfile)
  plot(bx)
}
}
\seealso{
\code{\link{bb_config}}, \code{\link{bb_example_sources}}, \code{\link{bb_source}}, \code{\link{bb_sync}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/config.R
\name{bb_config}
\alias{bb_config}
\title{Initialize a bowerbird configuration}
\usage{
bb_config(
  local_file_root,
  wget_global_flags = list(restrict_file_names = "windows", progress = "dot:giga"),
  http_proxy = NULL,
  ftp_proxy = NULL,
  clobber = 1
)
}
\arguments{
\item{local_file_root}{string: location of data repository on local file system}

\item{wget_global_flags}{list: wget flags that will be applied to all data sources that call \code{bb_wget}. These will be appended to the data-source-specific wget flags provided via the source's method argument}

\item{http_proxy}{string: URL of HTTP proxy to use e.g. 'http://your.proxy:8080' (NULL for no proxy)}

\item{ftp_proxy}{string: URL of FTP proxy to use e.g. 'http://your.proxy:21' (NULL for no proxy)}

\item{clobber}{numeric: 0=do not overwrite existing files, 1=overwrite if the remote file is newer than the local copy, 2=always overwrite existing files}
}
\value{
configuration object
}
\description{
The configuration object controls the behaviour of the bowerbird synchronization process, run via \code{bb_sync(my_config)}. The configuration object defines the data sources that will be synchronized, where the data files from those sources will be stored, and a range of options controlling how the synchronization process is conducted. The parameters provided here are repository-wide settings, and will affect all data sources that are subsequently added to the configuration.
}
\details{
Note that the \code{local_file_root} directory need not actually exist when the configuration object is created, but when \code{bb_sync} is run, either the directory must exist or \code{create_root=TRUE} must be passed (i.e. \code{bb_sync(...,create_root=TRUE)}).
}
\examples{
\dontrun{
  cf <- bb_config("/my/file/root") \%>\%
    bb_add(bb_example_sources())

  ## save to file
  saveRDS(cf,file="my_config.rds")
  ## load previously saved config
  cf <- readRDS(file="my_config.rds")
}

}
\seealso{
\code{\link{bb_source}}, \code{\link{bb_sync}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/re_exports.R
\docType{import}
\name{reexports}
\alias{reexports}
\alias{\%>\%}
\title{Objects exported from other packages}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{magrittr}{\code{\link[magrittr:pipe]{\%>\%}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/config.R
\name{bb_settings}
\alias{bb_settings}
\alias{bb_settings<-}
\title{Gets or sets a bowerbird configuration object's settings}
\usage{
bb_settings(config)

bb_settings(config) <- value
}
\arguments{
\item{config}{bb_config: a bowerbird configuration (as returned by \code{bb_config})}

\item{value}{list: new values to set}
}
\value{
named list
}
\description{
Gets or sets a bowerbird configuration object's settings. These are repository-wide settings that are applied to all data sources added to the configuration. Use this function to alter the settings of a configuration previously created using \code{bb_config}.
}
\details{
Note that an assignment along the lines of \code{bb_settings(cf) <- new_settings} replaces all of the settings in the configuration with the \code{new_settings}. The most common usage pattern is to read the existing settings, modify them as needed, and then rewrite the whole lot back into the configuration object (as per the examples here).
}
\examples{
cf <- bb_config(local_file_root="/your/data/directory")

## see current settings
bb_settings(cf)

## add an http proxy
sets <- bb_settings(cf)
sets$http_proxy <- "http://my.proxy"
bb_settings(cf) <- sets

## change the current local_file_root setting
sets <- bb_settings(cf)
sets$local_file_root <- "/new/location"
bb_settings(cf) <- sets

}
\seealso{
\code{\link{bb_config}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sync.R
\name{bb_sync}
\alias{bb_sync}
\title{Run a bowerbird data repository synchronization}
\usage{
bb_sync(
  config,
  create_root = FALSE,
  verbose = FALSE,
  catch_errors = TRUE,
  confirm_downloads_larger_than = 0.1,
  dry_run = FALSE
)
}
\arguments{
\item{config}{bb_config: configuration as returned by \code{\link{bb_config}}}

\item{create_root}{logical: should the data root directory be created if it does not exist? If this is \code{FALSE} (default) and the data root directory does not exist, an error will be generated}

\item{verbose}{logical: if \code{TRUE}, provide additional progress output}

\item{catch_errors}{logical: if \code{TRUE}, catch errors and continue the synchronization process. The sync process works through data sources sequentially, and so if \code{catch_errors} is \code{FALSE}, then an error during the synchronization of one data source will prevent all subsequent data sources from synchronizing}

\item{confirm_downloads_larger_than}{numeric or NULL: if non-negative, \code{bb_sync} will ask the user for confirmation to download any data source of size greater than this number (in GB). A value of zero will trigger confirmation on every data source. A negative or NULL value will not prompt for confirmation. Note that this only applies when R is being used interactively. The expected download size is taken from the \code{collection_size} parameter of the data source, and so its accuracy is dependent on the accuracy of the data source definition}

\item{dry_run}{logical: if \code{TRUE}, \code{bb_sync} will do a dry run of the synchronization process without actually downloading files}
}
\value{
a tibble with the \code{name}, \code{id}, \code{source_url}, sync success \code{status}, and \code{files} of each data source. Data sources that contain multiple source URLs will appear as multiple rows in the returned tibble, one per \code{source_url}. \code{files} is a tibble with columns \code{url} (the URL the file was downloaded from), \code{file} (the path to the file), and \code{note} (either "downloaded" for a file that was downloaded, "local copy" for a file that was not downloaded because there was already a local copy, or "decompressed" for files that were extracted from a downloaded (or already-locally-present) compressed file. \code{url} will be \code{NA} for "decompressed" files
}
\description{
This function takes a bowerbird configuration object and synchronizes each of the data sources defined within it. Data files will be downloaded if they are not present on the local machine, or if the configuration has been set to update local files.
}
\details{
Note that when \code{bb_sync} is run, the \code{local_file_root} directory must exist or \code{create_root=TRUE} must be specified (i.e. \code{bb_sync(...,create_root=TRUE)}). If \code{create_root=FALSE} and the directory does not exist, \code{bb_sync} will fail with an error.
}
\examples{
\dontrun{
  ## Choose a location to store files on the local file system.
  ## Normally this would be an explicit choice by the user, but here
  ## we just use a temporary directory for example purposes.

  td <- tempdir()
  cf <- bb_config(local_file_root = td)

  ## Bowerbird must then be told which data sources to synchronize.
  ## Let's use data from the Australian 2016 federal election, which is provided as one
  ## of the example data sources:

  my_source <- bb_example_sources("Australian Election 2016 House of Representatives data")

  ## Add this data source to the configuration:

  cf <- bb_add(cf, my_source)

  ## Once the configuration has been defined and the data source added to it,
  ## we can run the sync process.
  ## We set \code{verbose=TRUE} so that we see additional progress output:

  status <- bb_sync(cf, verbose = TRUE)

  ## The files in this data set have been stored in a data-source specific
  ## subdirectory of our local file root:

  status$files[[1]]

  ## We can run this at any later time and our repository will update if the source has changed:

  status2 <- bb_sync(cf, verbose = TRUE)
}

}
\seealso{
\code{\link{bb_config}}, \code{\link{bb_source}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/config.R
\name{bb_data_sources}
\alias{bb_data_sources}
\alias{bb_data_sources<-}
\title{Gets or sets a bowerbird configuration object's data sources}
\usage{
bb_data_sources(config)

bb_data_sources(config) <- value
}
\arguments{
\item{config}{bb_config: a bowerbird configuration (as returned by \code{\link{bb_config}})}

\item{value}{data.frame: new data sources to set (e.g. as returned by \code{\link{bb_example_sources}})}
}
\value{
a tibble with columns as specified by \code{\link{bb_source}}
}
\description{
Gets or sets the data sources contained in a bowerbird configuration object.
}
\details{
Note that an assignment along the lines of \code{bb_data_sources(cf) <- new_sources} replaces all of the sources in the configuration with the \code{new_sources}. If you wish to modify the existing sources then read them, modify as needed, and then rewrite the whole lot back into the configuration object.
}
\examples{
## create a configuration and add data sources
cf <- bb_config(local_file_root="/your/data/directory")
cf <- bb_add(cf,bb_example_sources())

## examine the sources contained in cf
bb_data_sources(cf)

## replace the sources with different ones
\dontrun{
bb_data_sources(cf) <- new_sources
}

}
\seealso{
\code{\link{bb_config}}, \code{\link{bb_source}}, \code{\link{bb_example_sources}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/example_sources.R
\name{bb_source_us_buildings}
\alias{bb_source_us_buildings}
\title{Example bowerbird data source: Microsoft US Buildings}
\usage{
bb_source_us_buildings(states)
}
\arguments{
\item{states}{character: (optional) one or more US state names for which to download data. If missing, data from all states will be downloaded. See the reference page for valid state names}
}
\value{
a tibble with columns as specified by \code{\link{bb_source}}
}
\description{
This function constructs a data source definition for the Microsoft US Buildings data set. This data set contains 124,885,597 computer generated building footprints in all 50 US states. NOTE: currently, the downloaded zip files will not be unzipped automatically. Work in progress.
}
\examples{
\dontrun{
## define a configuration and add this buildings data source to it
##  only including data for the District of Columbia and Hawaii
cf <- bb_config(tempdir()) \%>\%
  bb_add(bb_source_us_buildings(states = c("District of Columbia", "Hawaii")))

## synchronize (download) the data
bb_sync(cf)
}

}
\references{
\url{https://github.com/Microsoft/USBuildingFootprints}
}
\seealso{
\code{\link{bb_example_sources}}, \code{\link{bb_config}}, \code{\link{bb_handler_rget}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/postprocess.R
\name{bb_cleanup}
\alias{bb_cleanup}
\title{Postprocessing: remove unwanted files}
\usage{
bb_cleanup(
  pattern,
  recursive = FALSE,
  ignore_case = FALSE,
  all_files = FALSE,
  ...
)
}
\arguments{
\item{pattern}{string: regular expression, passed to \code{file.info}}

\item{recursive}{logical: should the cleanup recurse into subdirectories?}

\item{ignore_case}{logical: should pattern matching be case-insensitive?}

\item{all_files}{logical: should the cleanup include hidden files?}

\item{...}{: extra parameters passed automatically by \code{bb_sync}}
}
\value{
a list, with components \code{status} (TRUE on success) and \code{deleted_files} (character vector of paths of files that were deleted)
}
\description{
A function for removing unwanted files after downloading. This function is not intended to be called directly, but rather is specified as a \code{postprocess} option in \code{\link{bb_source}}.
}
\details{
This function can be used to remove unwanted files after a data source has been synchronized. The \code{pattern} specifies a regular expression that is passed to \code{file.info} to find matching files, which are then deleted. Note that only files in the data source's own directory (i.e. its subdirectory of the \code{local_file_root} specified in \code{bb_config}) are subject to deletion. But, beware! Some data sources may share directories, which can lead to unexpected file deletion. Be as specific as you can with the \code{pattern} parameter.
}
\examples{
\dontrun{
  ## remove .asc files after synchronization
  my_source <- bb_source(..., postprocess = list(list("bb_cleanup", pattern = "\\\\.asc$")))
}

}
\seealso{
\code{\link{bb_source}}, \code{\link{bb_config}}, \code{\link{bb_decompress}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bowerbird.R
\docType{package}
\name{bowerbird}
\alias{bowerbird}
\alias{bowerbird-package}
\title{\pkg{bowerbird}}
\description{
Often it's desirable to have local copies of third-party data sets. Fetching data on the fly from remote sources can be a great strategy, but for speed or other reasons it may be better to have local copies. This is particularly common in environmental and other sciences that deal with large data sets (e.g. satellite or global climate model products). Bowerbird is an R package for maintaining a local collection of data sets from a range of data providers.
}
\references{
\url{https://github.com/AustralianAntarcticDivision/bowerbird}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/config.R
\name{bb_data_source_dir}
\alias{bb_data_source_dir}
\title{Return the local directory of each data source in a configuration}
\usage{
bb_data_source_dir(config)
}
\arguments{
\item{config}{bb_config: configuration as returned by \code{\link{bb_config}}}
}
\value{
character vector of directories
}
\description{
Return the local directory of each data source in a configuration. Files from each data source are stored locally in the associated directory. Note that if a data source has multiple \code{source_url} values, this function might return multiple directory names (depending on whether those \code{source_url}s map to the same directory or not).
}
\examples{
cf <- bb_config("/my/file/root") \%>\%
  bb_add(bb_example_sources())
bb_data_source_dir(cf)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oceandata.R
\name{bb_handler_oceandata}
\alias{bb_handler_oceandata}
\title{Handler for Oceandata data sets}
\usage{
bb_handler_oceandata(search, dtype, sensor, ...)
}
\arguments{
\item{search}{string: (required) the search string to pass to the oceancolor file searcher (https://oceandata.sci.gsfc.nasa.gov/api/file_search)}

\item{dtype}{string: (optional) the data type (e.g. "L3m") to pass to the oceancolor file searcher. Valid options at the time of writing are aquarius, seawifs, aqua, terra, meris, octs, czcs, hico, viirs (for snpp), viirsj1, s3olci (for sentinel-3a), s3bolci (see https://oceancolor.gsfc.nasa.gov/data/download_methods/)}

\item{sensor}{string: (optional) the sensor (e.g. "seawifs") to pass to the oceancolor file searcher. Valid options at the time of writing are L0, L1, L2, L3b (for binned data), L3m (for mapped data), MET (for ancillary data), misc (for sundry products)}

\item{...}{: extra parameters passed automatically by \code{bb_sync}}
}
\value{
TRUE on success
}
\description{
This is a handler function to be used with data sets from NASA's Oceandata system. This function is not intended to be called directly, but rather is specified as a \code{method} option in \code{\link{bb_source}}.
}
\details{
Note that users will need an Earthdata login, see https://urs.earthdata.nasa.gov/. Users will also need to authorize the application 'OB.DAAC Data Access' (see 'My Applications' at https://urs.earthdata.nasa.gov/profile)

Oceandata uses standardized file naming conventions (see https://oceancolor.gsfc.nasa.gov/docs/format/), so once you know which products you want you can construct a suitable file name pattern to search for. For example, "S*L3m_MO_CHL_chlor_a_9km.nc" would match monthly level-3 mapped chlorophyll data from the SeaWiFS satellite at 9km resolution, in netcdf format. This pattern is passed as the \code{search} argument. Note that the \code{bb_handler_oceandata} does not take need `source_url` to be specified in the \code{bb_source} call.
}
\examples{

my_source <- bb_source(
  name="Oceandata SeaWiFS Level-3 mapped monthly 9km chl-a",
  id="SeaWiFS_L3m_MO_CHL_chlor_a_9km",
  description="Monthly remote-sensing chlorophyll-a from the SeaWiFS satellite at
    9km spatial resolution",
  doc_url="https://oceancolor.gsfc.nasa.gov/",
  citation="See https://oceancolor.gsfc.nasa.gov/citations",
  license="Please cite",
  method=list("bb_handler_oceandata",search="S*L3m_MO_CHL_chlor_a_9km.nc"),
  postprocess=NULL,
  collection_size=7.2,
  data_group="Ocean colour")

}
\references{
https://oceandata.sci.gsfc.nasa.gov/
}
