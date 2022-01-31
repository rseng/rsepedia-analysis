rinat: Access iNaturalist data with R
================
Edmund Hart, Stéphane Guillou

[![Build
Status](https://api.travis-ci.org/ropensci/rinat.png)](https://travis-ci.org/ropensci/rinat)
[![Build
status](https://ci.appveyor.com/api/projects/status/gv7s9um107bep4na/branch/master)](https://ci.appveyor.com/project/sckott/rinat/branch/master)
[![codecov.io](https://codecov.io/github/ropensci/rinat/coverage.svg?branch=master)](https://codecov.io/github/ropensci/rinat?branch=master)
[![](https://cranlogs.r-pkg.org/badges/rinat)](https://CRAN.R-project.org/package=rinat)

R wrapper for iNaturalist APIs for accessing the observations. The
detailed documentation of the API is available on the [iNaturalist
website](https://www.inaturalist.org/pages/api+reference) and is part of
our larger species occurrence searching packages
[SPOCC](https://github.com/ropensci/spocc).

## Installation

You can install the latest version available on CRAN with:

``` r
install.packages("rinat")
```

Alternatively, you can install the development version from Github with:

``` r
remotes::install_github("ropensci/rinat")
```

## Usage

### Get observations

#### Text search

You can search for observations by either common or scientific name. It
will search the entire iNaturalist database, so the search below will
return all entries that *mention* Monarch butterflies, not just Monarch
observations.

``` r
library(rinat)
monarchs <- get_inat_obs(query = "Monarch Butterfly")
unique(monarchs$scientific_name)
```

    ## [1] "Danaus plexippus" "Danaina"

> Note that `get_inat_obs()` will return 100 observations by default.
> This can be controlled with the `maxresults` argument.

Another use for a fuzzy search is searching for a habitat,
e.g. searching for all observations that might happen in a vernal pool.
We can then see all the taxon names found.

``` r
vp_obs <- get_inat_obs(query = "vernal pool")
# see the first few taxa
head(vp_obs$scientific_name)
```

    ## [1] "Micrathetis triplex"     "Sphaeropthalma unicolor"
    ## [3] "Lepidoptera"             "Lepidoptera"            
    ## [5] "Synchlora faseolaria"    "Lepidoptera"

#### Taxon search

To return only records of a specific species or taxonomic group, use the
`taxon_name` argument. For example, to return observations of anything
from the Nymphalidae family, and restricting the search to the year
2015:

``` r
nymphalidae <- get_inat_obs(taxon_name  = "Nymphalidae", year = 2015)
# how many unique taxa?
length(unique(nymphalidae$scientific_name))
```

    ## [1] 79

And to return only the Monarch butterfly observations that also mention
the term “chrysalis”:

``` r
monarch_chrysalis <- get_inat_obs(taxon_name = "Danaus plexippus", query = "chrysalis")
```

#### Bounding box search

You can also search within a bounding box by giving a simple set of
coordinates.

``` r
## Search by area
bounds <- c(38.44047, -125, 40.86652, -121.837)
deer <- get_inat_obs(query = "Mule Deer", bounds = bounds)
plot(deer$longitude, deer$latitude)
```

![](README_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

### Other functions

More functions are available, notably to access:

  - observations in a project with `get_inat_obs_project()`
  - details of a single observation with `get_inat_obs_id()`
  - observations from a single user with `get_inat_obs_user()`
  - taxa statistics with `get_inat_taxon_stats()`
  - user statistics with `get_inat_user_stats()`

More detailed examples are included in the vignette:

``` r
vignette("rinat-intro", package = "rinat")
```

#### Mapping

Basic maps can be created as well to quickly visualize search results.
Maps can either be plotted automatically with `plot = TRUE` (the
default), or simply return a ggplot2 object with `plot = FALSE`. This
works well with single species data, but more complicated plots are best
made from scratch.

``` r
library(ggplot2)

## Map 100 spotted salamanders
a_mac <- get_inat_obs(taxon_name = "Ambystoma maculatum")
salamander_map <- inat_map(a_mac, plot = FALSE)

### Now we can modify the returned map
salamander_map + borders("state") + theme_bw()
```

<img src="README_files/figure-gfm/unnamed-chunk-9-1.png" width="672" />

`inat_map()` is useful for quickly mapping data obtained with rinat.
Here is an example of customised map that does not make use of it. (Not
the use of `quality = "research"` to restrict the search to the more
reliable observations.)

``` r
## A more elaborate map of Colibri sp.
colibri <- get_inat_obs(taxon_name = "Colibri",
                        quality = "research",
                        maxresults = 500)
ggplot(data = colibri, aes(x = longitude,
                         y = latitude,
                         colour = scientific_name)) +
  geom_polygon(data = map_data("world"),
                   aes(x = long, y = lat, group = group),
                   fill = "grey95",
                   color = "gray40",
                   size = 0.1) +
  geom_point(size = 0.7, alpha = 0.5) +
  coord_fixed(xlim = range(colibri$longitude, na.rm = TRUE),
              ylim = range(colibri$latitude, na.rm = TRUE)) +
  theme_bw()
```

<img src="README_files/figure-gfm/unnamed-chunk-10-1.png" width="672" />

-----

[![](http://ropensci.org/public_images/github_footer.png)](https://ropensci.org/)
# rinat 0.1.8

* Properly cater for mismatch in project observation numbers to stay on CRAN (attempt in previous version did not cover all cases)
* Improve console messages in `get_inat_obs_project()`
* Skip tests that use the iNaturalist API on CRAN
* Don't error when iNaturalist or Internet is not available, as per CRAN Repository Policies (#49)

# rinat 0.1.7

* Cater for a corner case in which a mismatch between a project's reported number of observations and the actual number of observations returned produced an error when merging data frames in `get_inat_obs_project()`

# rinat 0.1.6

* New maintainer: Stéphane Guillou (@stragu)
* Improve documentation
* Clean up code according to CRAN feedback
* Note that 0.1.x versions of rinat will only fix existing issues, before a move to the new iNaturalist API with version 0.2 (which is likely to introduce breaking changes).

## New features

* `get_inat_obs()` can now use objects of class `Spatial` or `sf` as a bounding box (#38, @Martin-Jung and #40, @LDalby; fixes #30)
* Allow the use of an iNat `place_id` when searching for observations (commit 1d7f14f, @LDalby)

## Bug fixes

* Lower result number limit when bounding box present (#20)
* Fix result pagination in `get_obs_project()` for cases when the number of results is a multiple of 100 (#1, @vijaybarve)
* Stop `get_inat_obs()` if no search parameter(s) provided (#32, @LDalby)
* Avoid API error by limiting search results to 10000 (#32, @LDalby)
* Fix argument name in `inat_handle()` (#32, @LDalby)
* Use JSON endpoint instead of CSV to avoid `get_inat_obs_project()` failing on some projects (#35, @LDalby and #38, @Martin-Jung; fixes #28 and #37)
* Fix code according to `devtools::check()` feedback in preparation for CRAN release

# rinat 0.1.5

## Bug fixes

* Fixed bug where an error occurred when >20K records were requested and code now throws a warning to not hammer the API
* Fixed warning thrown when building vignettes on Linux systems
* Fixed bug where example code parameter names were different than actual parameter names

## New features

* Added NEWS file.
* Added a full suite of tests
* Added new vignette that builds with markdown and not hacky prebuilt PDF## Test environments

* local Ubuntu 18.04, R 4.0.4
* win-builder R-devel with `devtools::check_win_devel()`
* win-builder R-release with `devtools::check_win_release()`

## R CMD check results

There were no ERRORs, no WARNINGs

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'St�phane Guillou <stephane.guillou@member.fsf.org>'
  New submission
  Package was archived on CRAN
  Possibly mis-spelled words in DESCRIPTION:
  APIs (3:42)
  CRAN repository db overrides:
  X-CRAN-Comment: Archived on 2021-03-04 for policy violation.
  
Misspelling is irrelevant, and reason for archival is addressed in this release.---
title: "rinat: Access iNaturalist data with R"
author: Edmund Hart, Stéphane Guillou
output: github_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, warning = FALSE)
```

[![Build Status](https://api.travis-ci.org/ropensci/rinat.png)](https://travis-ci.org/ropensci/rinat)
[![Build status](https://ci.appveyor.com/api/projects/status/gv7s9um107bep4na/branch/master)](https://ci.appveyor.com/project/sckott/rinat/branch/master)
[![codecov.io](https://codecov.io/github/ropensci/rinat/coverage.svg?branch=master)](https://codecov.io/github/ropensci/rinat?branch=master)
[![](https://cranlogs.r-pkg.org/badges/rinat)](https://CRAN.R-project.org/package=rinat)


R wrapper for iNaturalist APIs for accessing the observations. The detailed documentation of the API is available on the [iNaturalist website](https://www.inaturalist.org/pages/api+reference) and is part of our larger species occurrence searching packages [SPOCC](https://github.com/ropensci/spocc).


## Installation

You can install the latest version available on CRAN with:

```{r eval=FALSE}
install.packages("rinat")
```

Alternatively, you can install the development version from Github with:

```{r eval=FALSE}
remotes::install_github("ropensci/rinat")
```

## Usage

### Get observations

#### Text search

You can search for observations by either common or scientific name. It will search the entire iNaturalist database, so the search below will return all entries that _mention_ Monarch butterflies, not just Monarch observations.

```{r}
library(rinat)
monarchs <- get_inat_obs(query = "Monarch Butterfly")
unique(monarchs$scientific_name)
```
> Note that `get_inat_obs()` will return 100 observations by default. This can be controlled with the `maxresults` argument.

Another use for a fuzzy search is searching for a habitat, e.g. searching for all observations that might happen in a vernal pool. We can then see all the taxon names found.

```{r}
vp_obs <- get_inat_obs(query = "vernal pool")
# see the first few taxa
head(vp_obs$scientific_name)
```


#### Taxon search

To return only records of a specific species or taxonomic group, use the `taxon_name` argument. For example, to return observations of anything from the Nymphalidae family, and restricting the search to the year 2015:

```{r}
nymphalidae <- get_inat_obs(taxon_name  = "Nymphalidae", year = 2015)
# how many unique taxa?
length(unique(nymphalidae$scientific_name))
```

And to return only the Monarch butterfly observations that also mention the term "chrysalis":

```{r}
monarch_chrysalis <- get_inat_obs(taxon_name = "Danaus plexippus", query = "chrysalis")
```

#### Bounding box search

You can also search within a bounding box by giving a simple set of coordinates.

```{r}
## Search by area
bounds <- c(38.44047, -125, 40.86652, -121.837)
deer <- get_inat_obs(query = "Mule Deer", bounds = bounds)
plot(deer$longitude, deer$latitude)
```


### Other functions

More functions are available, notably to access:

* observations in a project with `get_inat_obs_project()`
* details of a single observation with `get_inat_obs_id()`
* observations from a single user with `get_inat_obs_user()`
* taxa statistics with `get_inat_taxon_stats()`
* user statistics with `get_inat_user_stats()`

More detailed examples are included in the vignette:

```{r eval=FALSE}
vignette("rinat-intro", package = "rinat")
```


#### Mapping

Basic maps can be created as well to quickly visualize search results. Maps can either be plotted automatically with `plot = TRUE` (the default), or simply return a ggplot2 object with `plot = FALSE`. This works well with single species data, but more complicated plots are best made from scratch.

```{r fig.width=7, fig.height=4, fig.retina=3}
library(ggplot2)

## Map 100 spotted salamanders
a_mac <- get_inat_obs(taxon_name = "Ambystoma maculatum")
salamander_map <- inat_map(a_mac, plot = FALSE)

### Now we can modify the returned map
salamander_map + borders("state") + theme_bw()
```

`inat_map()` is useful for quickly mapping data obtained with rinat. Here is an example of customised map that does not make use of it. (Not the use of `quality = "research"` to restrict the search to the more reliable observations.)

```{r fig.width=7, fig.height=7, fig.retina=3}
## A more elaborate map of Colibri sp.
colibri <- get_inat_obs(taxon_name = "Colibri",
                        quality = "research",
                        maxresults = 500)
ggplot(data = colibri, aes(x = longitude,
                         y = latitude,
                         colour = scientific_name)) +
  geom_polygon(data = map_data("world"),
                   aes(x = long, y = lat, group = group),
                   fill = "grey95",
                   color = "gray40",
                   size = 0.1) +
  geom_point(size = 0.7, alpha = 0.5) +
  coord_fixed(xlim = range(colibri$longitude, na.rm = TRUE),
              ylim = range(colibri$latitude, na.rm = TRUE)) +
  theme_bw()
```


---

[![](http://ropensci.org/public_images/github_footer.png)](https://ropensci.org/)---
title: "Access iNaturalist data through APIs"
author: "Edmund Hart, Stéphane Guillou"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Access iNaturalist data through APIs}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, warning = FALSE)
```

## About

rinat is a wrapper for iNaturalist APIs for accessing the observations. The detailed documentation of the API is available on the [iNaturalist website](https://www.inaturalist.org/pages/api+reference) and is part of our larger species occurrence searching packages [SPOCC](https://github.com/ropensci/spocc).


## Quickstart guide

### Get observations

#### Fuzzy search

You can search for observations by either common or scientific name. It will search the entire iNaturalist database, so the search below will return all entries that mention Monarch butterflies, not just entries for Monarchs.

```{r}
library(rinat)
butterflies <- get_inat_obs(query = "Monarch Butterfly")
```

Another use for a fuzzy search is searching for a common name or habitat, e.g. searching for all observations that might happen in a vernal pool. We can then see all the species names found.

```{r}
vp_obs <- get_inat_obs(query = "vernal pool")
head(vp_obs$species_guess)
```


#### Taxon query

To return only records for a specific species or taxonomic group, use the taxon option.

```{r}
## Return observations in the family Nymphalidae, for 2015 only
nymphalidae <- get_inat_obs(taxon_name  = "Nymphalidae", year = 2015)

## Return just Monarch Butterfly records, all years
monarchs <- get_inat_obs(taxon_name = "Danaus plexippus")
```


#### Bounding box search

You can also search within a bounding box by giving a simple set of coordinates.

```{r fig.width=7, fig.height=4, fig.retina=3}
## Search by area
bounds <- c(38.44047, -125, 40.86652, -121.837)
deer <- get_inat_obs(query = "Mule Deer", bounds = bounds)
plot(deer$longitude, deer$latitude)
```


### Other functions

#### Get information and observations by project

You can get all the observations for a project if you know its ID or name as an iNaturalist slug.

```{r}
## Just get info about a project
vt_crows <- get_inat_obs_project("crows-in-vermont", type = "info", raw = FALSE)
```

```{r}
## Now get all the observations for that project
vt_crows_obs <- get_inat_obs_project(vt_crows$id, type = "observations")
```


#### Get observation details

Detailed information about a specific observation can be retrieved by observation ID. The easiest way to get the ID is from a more general search.

```{r}
m_obs <- get_inat_obs(query = "Monarch Butterfly")
head(get_inat_obs_id(m_obs$id[1]))
```


#### Get all observations by user

If you just want all the observations by a user you can download all their observations by user ID. A word of warning though, this can be quite large (easily into the 1000's).

```{r}
user_obs <- get_inat_obs_user(m_obs$user_login[1], maxresults = 20)
head(user_obs)[,1:5]
```


#### Stats by taxa

Basic statistics are available for taxa counts by date, date range, place ID (numeric ID), or user ID (string). Only the top 5 species are listed.

```{r}
## By date
counts <- get_inat_taxon_stats(date = "2020-06-14")
counts$total
### Top 5 species
counts$species_counts
### Most common taxon ranks
counts$rank_counts
```


#### Stats by user

Similar statistics can be gotten for users. The same input parameters can be used.

```{r}
## By date
counts <- get_inat_user_stats(date = "2010-06-14")
counts$total
counts$most_observations[1:10,]
counts$most_species[1:10,]
```

```{r}
## By place_ID
vt_crows <- get_inat_obs_project("crows-in-vermont", type = "info", raw = FALSE)
place_counts <- get_inat_user_stats(place = vt_crows$place_id)
place_counts$total
place_counts$most_observations[1:10,]
place_counts$most_species[1:10,]
```


### Mapping

Basic maps can be created as well to quickly visualize search results. Maps can either be plotted automatically with `plot = TRUE` (the default), or simply return a ggplot2 object with `plot = FALSE`.  This works well with single species data, but more complicated plots are best made from scratch.

```{r fig.width=7, fig.height=4, fig.retina=3}
library(ggplot2)

## Map 100 spotted salamanders
a_mac <- get_inat_obs(taxon_name = "Ambystoma maculatum")
salamander_map <- inat_map(a_mac, plot = FALSE)

### Now we can modify the returned map
salamander_map + borders("state") + theme_bw()
```

```{r fig.width=7, fig.height=7, fig.retina=3}
## A more elaborate map of Colibri sp.
colibri <- get_inat_obs(taxon_name = "Colibri",
                        quality = "research",
                        maxresults = 500)
ggplot(data = colibri, aes(x = longitude,
                         y = latitude,
                         colour = scientific_name)) +
  geom_polygon(data = map_data("world"),
                   aes(x = long, y = lat, group = group),
                   fill = "grey95",
                   color = "gray40",
                   size = 0.1) +
  geom_point(size = 0.7, alpha = 0.5) +
  coord_fixed(xlim = range(colibri$longitude, na.rm = TRUE),
              ylim = range(colibri$latitude, na.rm = TRUE)) +
  theme_bw()
```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_inat_user_stats.R
\name{get_inat_user_stats}
\alias{get_inat_user_stats}
\title{Get stats on users}
\usage{
get_inat_user_stats(
  date = NULL,
  date_range = NULL,
  place = NULL,
  project = NULL,
  uid = NULL
)
}
\arguments{
\item{date}{retrieve observations on a specific date, must be a string in the form YYYY-MM-DD}

\item{date_range}{a vector of dates, in the form YYYY-MM-DD}

\item{place}{get taxon stats by place, you can find place id's on the iNaturalist page: https://www.inaturalist.org/places, must be a numeric ID}

\item{project}{get taxon stats by project id}

\item{uid}{get taxon stats by user id}
}
\value{
a list with two data frames with of the 5 users with the most observations and the most species
}
\description{
Get stats on which users reported the most species or had tho most observations within a given range. This range can be by user, place, project, day or date range. Output will be a count of the total number of taxa observed at each taxonomic level.
}
\examples{
\dontrun{
 counts <- get_inat_user_stats(date = "2010-06-14")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_inat_obs_project.R
\name{get_inat_obs_project}
\alias{get_inat_obs_project}
\title{Download observations or info from a project}
\usage{
get_inat_obs_project(grpid, type = c("observations", "info"), raw = FALSE)
}
\arguments{
\item{grpid}{Name of the group as an iNaturalist slug or group id}

\item{type}{character Either "observations" or "info"  Observations returns all observations, and "info" returns project details similar to what you can find on a project webpage.}

\item{raw}{logical TRUE or FALSE. If TRUE and searching for project info, returns the raw output of parsed JSON for that project. Otherwise just some basic information is returned as a list}
}
\description{
retrieve observations from a particular iNaturalist project. This function can be used to get either observations or information from a project by project name or ID
}
\details{
An iNaturalist slug is usually the project as single string with words separated by hyphens. For instance, the project "State Flowers of the United States" has a slug of "state-flowers-of-the-united-states-eol-collection".  This can be extracted from the URL for the project usually. The state flowers project has the following URL http://www.inaturalist.org/projects/state-flowers-of-the-united-states-eol-collection
}
\examples{
\dontrun{
 get_inat_obs_project(354, type = "observations")
 get_inat_obs_project("crows-in-vermont", type="info",raw=FALSE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_inat_obs.R
\name{get_inat_obs}
\alias{get_inat_obs}
\title{Download iNaturalist data}
\usage{
get_inat_obs(
  query = NULL,
  taxon_name = NULL,
  taxon_id = NULL,
  place_id = NULL,
  quality = NULL,
  geo = NULL,
  year = NULL,
  month = NULL,
  day = NULL,
  bounds = NULL,
  maxresults = 100,
  meta = FALSE
)
}
\arguments{
\item{query}{Query string for a general search.}

\item{taxon_name}{Filter by iNat taxon name. Note that this will also select observations of
descendant taxa. Note that names are not unique, so if the name matches multiple taxa, no
observations may be returned.}

\item{taxon_id}{Filter by iNat taxon ID. Note that this will also select observations of descendant taxa.}

\item{place_id}{Filter by iNat place ID.}

\item{quality}{The quality grade to be used.  Must be either "casual" or "research".  If left
blank both will be returned.}

\item{geo}{Flag for returning only results that are georeferenced, TRUE will exclude
non-georeferenced results, but they cannot be excluded.}

\item{year}{Return observations only in that year (can only be one year, not a range of years).}

\item{month}{Return observations only by month, must be numeric, 1...12}

\item{day}{Return observations only on a given day of the month,  1...31}

\item{bounds}{A bounding box of longitude (-180 to 180) and latitude (-90 to 90) to search
within.  It is a vector in the form of southern latitude, western longitude, northern latitude,
and eastern longitude. Alternatively supply an sf or sp object from which the bounding box will
be derived.}

\item{maxresults}{The maximum number of results to return. Should not be
a number higher than 10000.}

\item{meta}{(logical) If TRUE, the output of this function is a list with metadata on the output
and a data.frame of the data. If FALSE (default), just the data.frame.}
}
\value{
A dataframe of the number of observations requested.
}
\description{
Primary function to retrieve observations from iNaturalist, allows users to search
for data, or just filter results by a subset of what is offered by the API.
}
\note{
Filtering doesn't always work with the query parameter for some reason (a problem on
the API end).  If you want to filter by time, it's best to use the scientific name and put it
in the 'taxa' field, and not in the query field.  Another issue is that the query parameter
will search the entire entry, so it is possible to get unintended results.  Depending on your
use case it may be advisable to use the "taxon" field instead of the query field.
}
\examples{
\dontrun{
  ### Make a standard query
  get_inat_obs(query = "Monarch Butterfly")

  ##Filter by a bounding box of Northern California
  bounds <- c(38.44047, -125, 40.86652, -121.837)
  get_inat_obs(query = "Mule Deer", bounds = bounds)

  ## Filter with by just taxon, allows higher order filtering,
  ## Here we can search for just stone flies (order Plecoptera)
  get_inat_obs(taxon_name = "Plecoptera")

  ## get metadata (the number of results found on the server)
  out <- get_inat_obs(query = "Monarch Butterfly", meta = TRUE)
  out$meta
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_inat_obs_user.R
\name{get_inat_obs_user}
\alias{get_inat_obs_user}
\title{Download observations for a user}
\usage{
get_inat_obs_user(username, maxresults = 100)
}
\arguments{
\item{username}{username of the iNaturalist user to fetch records}

\item{maxresults}{the maximum number of results to return}
}
\value{
a list with full details on a given record
}
\description{
Get all the observations of a specific iNaturalist user.
}
\examples{
\dontrun{
  m_obs <- get_inat_obs(query="Monarch Butterfly")
  get_inat_obs_user(as.character(m_obs$user_login[1]))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inat_map.R
\name{inat_map}
\alias{inat_map}
\title{Plot iNaturalist observations}
\usage{
inat_map(data, map = "usa", subregion = ".", plot = TRUE)
}
\arguments{
\item{data}{data frame of iNaturalist observations.}

\item{map}{the map region to plot, you can find full documentation in the \code{\link{map}} package, default is USA.}

\item{subregion}{the name of the subregion to plot, see full documentation in the \code{\link{map}} package.}

\item{plot}{a logical value. TRUE plots the map object and returns it, and FALSE returns a ggplot object that you can modify and plot later.}
}
\value{
A ggplot map object.
}
\description{
Plot observations from iNaturalist. You have the option of automatically plotting, or returning a ggplot map object that you can add layers onto.
}
\examples{
\dontrun{
  m_obs <- get_inat_obs(taxon_name = "Ambystoma maculatum")
  salamander_map <- inat_map(m_obs, plot = FALSE)
  ### Now we can modify the returned map
  salamander_map + borders("state") + theme_bw()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_inat_taxon_stats.R
\name{get_inat_taxon_stats}
\alias{get_inat_taxon_stats}
\title{Get stats on taxon counts}
\usage{
get_inat_taxon_stats(
  date = NULL,
  date_range = NULL,
  place = NULL,
  project = NULL,
  uid = NULL
)
}
\arguments{
\item{date}{retrieve observations on a specific date, must be a string in the form YYYY-MM-DD}

\item{date_range}{a vector of dates, in the form YYYY-MM-DD}

\item{place}{get taxon stats by place, you can find place id's on the iNaturalist page: https://www.inaturalist.org/places, must be a numeric ID}

\item{project}{get taxon stats by project id (numeric ID)}

\item{uid}{get taxon stats by user id (string)}
}
\value{
a vector listing counts of observations at each level of identification possible (species, genus, etc.)
}
\description{
Get stats on taxa within a constrained range. This range can be by user, place, project, day or date range. Output will be a count of the total number of taxa observed at each taxonomic level.
}
\examples{
\dontrun{
 counts <- get_inat_taxon_stats(date = "2010-06-14")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_inat_obs_id.R
\name{get_inat_obs_id}
\alias{get_inat_obs_id}
\title{Get information on a specific observation}
\usage{
get_inat_obs_id(id)
}
\arguments{
\item{id}{a single id for an iNaturalist observation record}
}
\value{
a list with full details on a given record
}
\description{
Get information on a specific observation
}
\examples{
\dontrun{
m_obs <- get_inat_obs(query="Monarch Butterfly")
get_inat_obs_id(m_obs$id[1])
}
}
