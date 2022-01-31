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

PostcodesioR
================

# PostcodesioR <img src='man/figures/logo.png' align="right" height="139" />

[![Travis-CI Build
Status](https://travis-ci.org/ropensci/PostcodesioR.svg?branch=master)](https://travis-ci.org/ropensci/PostcodesioR)
[![Package-License](https://img.shields.io/badge/license-GPL--3-brightgreen.svg?style=flat)](https://www.gnu.org/licenses/gpl-3.0.html)
[![](https://badges.ropensci.org/176_status.svg)](https://github.com/ropensci/software-review/issues/176)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/PostcodesioR)](https://cran.r-project.org/package=PostcodesioR)
[![DOI](https://zenodo.org/badge/64221541.svg)](https://zenodo.org/badge/latestdoi/64221541)
![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/PostcodesioR)

An API wrapper around [postcodes.io](https://postcodes.io/) - free UK
postcode lookup and geocoder. This package helps to find and transform
information about UK administrative geography like postcodes, LSOA,
MSOA, constituencies, counties, wards, districts, CCG or NUTS.

The package is based exclusively on open data provided by postcodes.io.
PostcodesioR can be used by data scientists or social scientists working
with geocoded UK data. A common task when working with such data is
aggregating geocoded data on different administrative levels,
e.g. turning postcode-level data into counties or regions. This package
can help in achieving this and in many other cases when changing the
aggregation of geographic data is required.

## Installation

This package can be installed from GitHub (developmental version) or
CRAN (stable).

In order to install PostcodesioR use one of the following commands:

``` r
# stable version
install.packages("PostcodesioR")
```

or

``` r
# developmental version
if(!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("ropensci/PostcodesioR")
```

## Loading

Load the package by typing

``` r
library(PostcodesioR)
```

## Examples

Where possible, I tried to return a data frame. Unfortunately, a lot of
API calls return more complex data and in those cases it is safer to use
lists. The API limits the number of returned calls. Check functions’
documentation for more details.

For additional information about the returned data and the function
calls see the original [documentation](https://postcodes.io/docs).

The main function of this package provides information related to a
given postcode

``` r
lookup_result <- postcode_lookup("EC1Y8LX")

#overview
str(lookup_result)
```

    ## 'data.frame':    1 obs. of  35 variables:
    ##  $ postcode                       : chr "EC1Y 8LX"
    ##  $ quality                        : int 1
    ##  $ eastings                       : int 532544
    ##  $ northings                      : int 182128
    ##  $ country                        : chr "England"
    ##  $ nhs_ha                         : chr "London"
    ##  $ longitude                      : num -0.0909
    ##  $ latitude                       : num 51.5
    ##  $ european_electoral_region      : chr "London"
    ##  $ primary_care_trust             : chr "Islington"
    ##  $ region                         : chr "London"
    ##  $ lsoa                           : chr "Islington 023D"
    ##  $ msoa                           : chr "Islington 023"
    ##  $ incode                         : chr "8LX"
    ##  $ outcode                        : chr "EC1Y"
    ##  $ parliamentary_constituency     : chr "Islington South and Finsbury"
    ##  $ admin_district                 : chr "Islington"
    ##  $ parish                         : chr "Islington, unparished area"
    ##  $ admin_county                   : logi NA
    ##  $ admin_ward                     : chr "Bunhill"
    ##  $ ced                            : logi NA
    ##  $ ccg                            : chr "NHS North Central London"
    ##  $ nuts                           : chr "Haringey and Islington"
    ##  $ admin_district_code            : chr "E09000019"
    ##  $ admin_county_code              : chr "E99999999"
    ##  $ admin_ward_code                : chr "E05000367"
    ##  $ parish_code                    : chr "E43000209"
    ##  $ parliamentary_constituency_code: chr "E14000764"
    ##  $ ccg_code                       : chr "E38000240"
    ##  $ ccg_id_code                    : chr "93C"
    ##  $ ced_code                       : chr "E99999999"
    ##  $ nuts_code                      : chr "TLI43"
    ##  $ lsoa_code                      : chr "E01002704"
    ##  $ msoa_code                      : chr "E02000576"
    ##  $ lau2_code                      : chr "E09000019"

Check the
[vignette](https://docs.ropensci.org/PostcodesioR/articles/Introduction.html)
to see all functions in action.

## Notes

Currently, there is a limit to the number of API calls that can be made.
However, [postcodes.io](https://postcodes.io/) provides full list of
geolocation data that can be used locally without limitations. The
original data is sourced from [Office for National Statistics Data
Portal](https://geoportal.statistics.gov.uk/). That
[file](https://github.com/ideal-postcodes/postcodes.io/blob/master/latest)
is rather large so I didn’t include it in the package.

Go to the package’s [website](https://docs.ropensci.org/PostcodesioR/)
or to my [blog](https://walczak.org/tag/postcodesior/) for more
examples.

## Code of Conduct

Please note that this project is released with a [Contributor Code of
Conduct](https://github.com/ropensci/PostcodesioR/blob/master/CONDUCT.md).
By participating in this project you agree to abide by its terms.

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# PostcodesioR 0.3.1

# PostcodesioR 0.3.0

* `scottish_postcode_lookup` added.
* New fields (codes) added.
* README updated (hex logo and downloads).
* New tests.

# PostcodesioR 0.2.0

* `bulk_postcode_lookup` bug fixed.

# PostcodesioR 0.1.1

* Added a `NEWS.md` file to track changes to the package.
## Update (v 0.3.1)

* Minor bus fixes.
* Updated tests.
* Updated documentation.

## Test environments
* local ubuntu 20.04, R 4.0.3
* win-builder (devel and release). No errors, warnings or notes
* R-hub using `rhub::check_for_cran()`. One note related to the change of maintainer's email.
* R-hub tested on:
- Windows Server 2008 R2 SP1, R-devel, 32/64 bit,
- Ubuntu Linux 16.04 LTS, R-release, GCC
- Fedora Linux, R-devel, clang, gfortran

## Update (v 0.3)

* Minor bus fixes.
* `scottish_postcode_lookup` added
* New fields (codes) added
* README updated (hex logo and downloads)
* New tests
* URLs in README updated

## Test environments
* local ubuntu 20.04, R 4.0.3
* win-builder (devel and release). No errors, warnings or notes
* R-hub using `rhub::check_for_cran()`. No errors, warnings or notes
* R-hub tested on:
- Windows Server 2008 R2 SP1, R-devel, 32/64 bit,
- Ubuntu Linux 16.04 LTS, R-release, GCC
- Fedora Linux, R-devel, clang, gfortran

## R CMD check results

0 errors | 0 warnings | 0 notes

## Resubmission

This is a resubmission. In this version I have:

* Wrapped Postcodes.io in DESCRIPTION Title in quotation marks

* Shortened the title

* Replaced `dontrun` with `donttest` in examples

## Test environments
* local ubuntu 16.04, R 3.6.1
* ubuntu 14.04 (on travis-ci), R 3.6.0
* win-builder (devel and release)
* R-hub using `rhub::check_for_cran()`

## R CMD check results

0 errors | 0 warnings | 0 notes

## `rhub::check_for_cran()` results

* Error

Installation failed with PREPERROR on R-hub's Fedora Linux, R-devel, clang, gfortran (and other platforms), on account of failures to install required dependencies. This seems to be an issue outside of my control.

* Note

Windows Server 2008 R2 SP1, R-devel, 32/64 bit and Ubuntu Linux 16.04 LTS, R-release, GCC return a note:

The note regarded potential mis-spelled words. These words were spelled correctly.

Possibly mis-spelled words in DESCRIPTION:
  geocoding (16:22, 18:20)
  io (3:37)

* This is a new release.
---
title: "PostcodesioR"
output: rmarkdown::github_document
---

# PostcodesioR <img src='man/figures/logo.png' align="right" height="139" />

[![Travis-CI Build Status](https://travis-ci.org/ropensci/PostcodesioR.svg?branch=master)](https://travis-ci.org/ropensci/PostcodesioR)
[![Package-License](https://img.shields.io/badge/license-GPL--3-brightgreen.svg?style=flat)](https://www.gnu.org/licenses/gpl-3.0.html)
[![](https://badges.ropensci.org/176_status.svg)](https://github.com/ropensci/software-review/issues/176)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/PostcodesioR)](https://cran.r-project.org/package=PostcodesioR)
[![DOI](https://zenodo.org/badge/64221541.svg)](https://zenodo.org/badge/latestdoi/64221541)
![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/PostcodesioR)

An API wrapper around [postcodes.io](https://postcodes.io/) - free UK postcode lookup and geocoder. This package helps to find and transform information about UK administrative geography like postcodes, LSOA, MSOA, constituencies, counties, wards, districts, CCG or NUTS.

The package is based exclusively on open data provided by postcodes.io. PostcodesioR can be used by data scientists or social scientists working with geocoded UK data. A common task when working with such data is aggregating geocoded data on different administrative levels, e.g. turning postcode-level data into counties or regions. This package can help in achieving this and in many other cases when changing the aggregation of geographic data is required.

## Installation

This package can be installed from GitHub (developmental version) or CRAN (stable).

In order to install PostcodesioR use one of the following commands:

```{r, eval = FALSE}
# stable version
install.packages("PostcodesioR")
```

or

```{r, eval = FALSE}
# developmental version
if(!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("ropensci/PostcodesioR")
```

## Loading

Load the package by typing

```{r, warning = FALSE, message = FALSE}
library(PostcodesioR)
```

## Examples

Where possible, I tried to return a data frame. Unfortunately, a lot of API calls return more complex data and in those cases it is safer to use lists. The API limits the number of returned calls. Check functions' documentation for more details.

For additional information about the returned data and the function calls see the original [documentation](https://postcodes.io/docs).

The main function of this package provides information related to a given postcode

```{r, message = FALSE, warning = FALSE}
lookup_result <- postcode_lookup("EC1Y8LX")

#overview
str(lookup_result)
```

Check the [vignette](https://docs.ropensci.org/PostcodesioR/articles/Introduction.html) to see all functions in action.

## Notes

Currently, there is a limit to the number of API calls that can be made. However, [postcodes.io](https://postcodes.io/) provides full list of geolocation data that can be used locally without limitations. The original data is sourced from [Office for National Statistics Data Portal](https://geoportal.statistics.gov.uk/).
That [file](https://github.com/ideal-postcodes/postcodes.io/blob/master/latest) is rather large so I didn't include it in the package.

Go to the package's [website](https://docs.ropensci.org/PostcodesioR/) or to my [blog](https://walczak.org/tag/postcodesior/) for more examples.

## Code of Conduct

Please note that this project is released with a [Contributor Code of Conduct](https://github.com/ropensci/PostcodesioR/blob/master/CONDUCT.md). By participating in this project you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
---
title: "Introduction"
author: "Eryk Walczak"
date: "`r Sys.Date()`"
output: 
  rmarkdown::html_vignette:
    smart: no
    toc: true
vignette: >
  %\VignetteIndexEntry{Introduction}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

PostcodesioR is an API wrapper for postcodes.io. It allows acquiring geographic information about the UK postcodes and geographic coordinates.

## Installation

```{r, message = FALSE, warning = FALSE, eval = FALSE}
if (!require("devtools")) install.packages("devtools")
devtools::install_github("ropensci/PostcodesioR")
```


## Lookup postcodes and outcodes

### Single postcode

Provide a postcode to obtain all available information 

```{r, message = FALSE, warning = FALSE}
library(PostcodesioR)

lookup_result <- postcode_lookup("EC1Y8LX")

#overview
str(lookup_result)
```

There is another function that returns the same data points but returns a list and allows optional parameters

```{r}
query_result <- postcode_query("EC1Y8LX")

#overview
str(query_result)
```

This function creates a nested list with the codes for administrative district, county, ward, parish, parliamentary constituency, CCG, and NUTS.

### Multiple postcodes

To query two or more postcodes, use `bulk_` functions.

```{r}
pc_list <- list(postcodes = c("PR3 0SG", "M45 6GN", "EX165BL"))
bulk_lookup_result <- bulk_postcode_lookup(pc_list)

#overview
str(bulk_lookup_result[1])
```

If you want to work with data frame then the nested list created above can be turned into a data frame

```{r}
library(purrr)

bulk_list <- lapply(bulk_lookup_result, "[[", 2)

bulk_df <-
  map_dfr(bulk_list,
          `[`,
          c("postcode", "longitude", "latitude"))
```

Querying Scottish postcodes requires a separate function:

```{r}
scottish_lookup <- scottish_postcode_lookup("EH12NG")

str(scottish_lookup)
```


### Outward code lookup

Provide an outcode to obtain geolocation data for the centroid of the specified outcode:

```{r}
ocl <- outward_code_lookup("E1")

#overview
str(ocl)
```

## Reverse geocoding

Provide latitude and longitude to obtain geographic information. Different levels of aggregation are available, i.e. postcode or outcode.

### Single postcode

```{r}
rev_geo <- reverse_geocoding(0.127, 51.507)

# overview
str(rev_geo[1])
```

### Multiple postcodes

To reverse geocode multiple values use the function underneath. The result is a nested list, which might be a bit intimidating, but it allows storing unequal number of elements.

```{r}
# create a list with the coordinates
geolocations_list <- structure(
 list(
 geolocations = structure(
 list(
 longitude = c(-3.15807731271522, -1.12935802905177),
 latitude = c(51.4799900627036, 50.7186356978817),
 limit = c(NA, 100L),
 radius = c(NA, 500L)),
 .Names = c("longitude", "latitude", "limit", "radius"),
 class = "data.frame",
 row.names = 1:2)),
 .Names = "geolocations")

bulk_rev_geo <- bulk_reverse_geocoding(geolocations_list)

bulk_rev_geo[[1]]$result[[1]]
```

The list above is not the most common way of storing files. It's more likely that a data frame will be used to store the geodata. In that case, it has to be turned into a list of a specific format required by the API:

```{r}
geolocations_df <- structure(
  list(
    longitude = c(-3.15807731271522, -1.12935802905177),
    latitude = c(51.4799900627036, 50.7186356978817),
    limit = c(NA, 100L),
    radius = c(NA, 500L)),
  .Names = c("longitude", "latitude", "limit", "radius"),
  row.names = 1:2,
  class = "data.frame")

geolocations_df

# turn a data frame into a list
geolocations_df2list <- list(geolocations_df)

# add a list name
names(geolocations_df2list) <- "geolocations"

# display correct input for the function
geolocations_df2list
```


Common usage of this function might be extracting particular variables. You can extract one variable like this:

```{r}
# extract one postcode
bulk_rev_geo[[1]]$result[[8]]$postcode
```

But more likely you will want more than one result. After all, that's the point of using a bulk function:

```{r}
# function to extract variables of interest
extract_bulk_geo_variable <- function(x) {
  bulk_results <- lapply(bulk_rev_geo, `[[`, "result")
  sapply(unlist(bulk_results, recursive = FALSE), `[[`, x)
}

# define the variables you need
variables_of_interest <- c("postcode", "latitude", "longitude")

# return a data frame with the variables
data.frame(
  sapply(variables_of_interest, extract_bulk_geo_variable))
```


### Single outcode

```{r}
out_rev_geocode <- outcode_reverse_geocoding("-3.15", "51.47")
# overview
str(out_rev_geocode[1])
```

## Generate random entries

### Postcodes

Generates a list with a random UK postcode and corresponding geographic information:

```{r}
# without restrictions
random_postcode()
```

A randomly generated postcode can also belong to a particular outcode:

```{r}
# restrict to an outcode
random_postcode("N1")
```

## Places

You can also generate a random place, specified by an OSGB code, with corresponding geographic information:

```{r}
random_place()
```


## Postcode validation

This function can validate a UK postcode:

```{r}
postcode_validation("EC1Y8LX") # actual UK postcode
```


```{r}
postcode_validation("XYZ") # incorrect UK postcode
```


## Autocomplete postcodes

Find the potential candidates for a postcode if you only know the beginning characters

```{r}
postcode_autocomplete("EC1")
```

It defaults to 10 candidates, but can be changed by specifying the `limit` argument. 

## Find nearest postcodes or outcodes

Provide a postcode to get a list of the nearest postcodes:

```{r}
near_pc <- nearest_postcode("EC1Y8LX")

#overview
str(near_pc[1])
```

You can also use outcodes:

```{r}
near_outcode <- nearest_outcode("EC1Y")

# overview
str(near_outcode[2])
```

Or longitude and latitude

```{r}
near_ll <- nearest_outcode_lonlat(0.127, 51.507)

#overview
str(near_ll[1])
```


## Find places

Provide a name of a place of interest. You can specify the number of results (default is 10):

```{r}
place_query_result <- place_query("Hills", limit = 11)

# overview
str(place_query_result[1])
```

You can also find a place using an OSGB code:

```{r}
place_lookup_result <- place_lookup("osgb4000000074544700")

# overview
str(place_lookup_result)
```

## Terminated postcodes

You might end up having terminated postcodes in your data set. These are postcodes that are no longer active. UK postcodes can change so it's worth checking whether used postcodes are still active. If you need more information about when a particular postcode was terminated use:

```{r}
terminated_postcode("E1W 1UU")
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nearest_postcode.R
\name{nearest_postcode}
\alias{nearest_postcode}
\title{Nearest postcode}
\usage{
nearest_postcode(postcode, limit = 10, radius = 100)
}
\arguments{
\item{postcode}{A string. Valid UK postcode.}

\item{limit}{A string or integer. Limits number of postcodes matches to return. Defaults to 10. Needs to be lower than 100.}

\item{radius}{Limits number of postcodes matches to return. Defaults to 100m. Needs to be less than 2,000m.}
}
\value{
A list of geographic properties of the nearest postcode.
}
\description{
Returns nearest postcodes for a given postcode.
The search is based on the relative distance of the postcode centroid.
}
\examples{
\donttest{
nearest_postcode("EC1Y 8LX")
nearest_postcode("EC1Y 8LX", limit = 11)
nearest_postcode("EC1Y 8LX", limit = 12, radius = 200)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/terminated_postcode.R
\name{terminated_postcode}
\alias{terminated_postcode}
\title{Terminated postcode lookup}
\usage{
terminated_postcode(postcode)
}
\arguments{
\item{postcode}{A string. Terminated UK postcode.}
}
\value{
A data frame with data about terminated postcode. NULL if postcode is active.
\itemize{
\item \code{postcode} Postcode. All currently terminated postcodes within
the United Kingdom, the Channel Islands and the Isle of Man,
received every 3 months from Royal Mail. 2, 3 or 4-character outward code,
single space and 3-character inward code.
\item \code{year_terminated} Termination year.
Year of termination of a postcode.
\item \code{month_terminated} Termination month.
Month of termination of a postcode. 1-January, 2-February, ..., 12-December.
\item \code{longitude} Longitude.
The WGS84 longitude given the Postcode's national grid reference.
\item \code{latitude} Latitude.
The WGS84 latitude given the Postcode's national grid reference.
}

See
\url{https://postcodes.io/docs} for more details.
}
\description{
Returns month and year if a postcode was terminated or is no longer active.
}
\examples{
\donttest{
terminated_postcode("EC1Y 8LX") # existing postcode
terminated_postcode("E1W 1UU") # terminated postcode
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scottish_postcode_lookup.R
\name{scottish_postcode_lookup}
\alias{scottish_postcode_lookup}
\title{Scottish postcode lookup}
\usage{
scottish_postcode_lookup(postcode)
}
\arguments{
\item{postcode}{A string. One valid Scottish postcode.
This function is case- and space-insensitive.
For non-Scottish postcodes use \code{\link{postcode_lookup}}.
For more than one non-Scottish postcode use \code{\link{bulk_postcode_lookup}}.}
}
\value{
A data frame. Returns all available data if found. Returns NAs if postcode does not exist (404).
\itemize{
\item \code{postcode} Postcode. Royal Mail postcode.
\item \code{scottish_parliamentary_constituency} Scottish Parliamentary
Constituency 2014 Scottish Parliamentary Constituency.
\item \code{scottish_parliamentary_constituency} Scottish Parliamentary
Constituency GSS Code. A code that identifies a 2014 Scottish
Parliamentary Constituency.
}

See
\url{https://postcodes.io/docs} for more details.
}
\description{
Lookup a Scottish postcode.
}
\examples{
\donttest{
scottish_postcode_lookup("EH12NG")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outcode_reverse_geocoding.R
\name{outcode_reverse_geocoding}
\alias{outcode_reverse_geocoding}
\title{Outcode reverse geocoding}
\usage{
outcode_reverse_geocoding(longitude, latitude, limit = 10, radius = 5000)
}
\arguments{
\item{longitude}{A string, integer or float. Needs to have at least two decimal points.}

\item{latitude}{A string, integer or float. Needs to have at least two decimal points.}

\item{limit}{A string, integer or float. Limits number of postcodes matches to return. Defaults to 10. Needs to be less than 100.}

\item{radius}{A string, integer or float. Limits number of postcodes matches to return. Defaults to 5,000m. Needs to be less than 25,000m.}
}
\value{
A list of geographical properties.
}
\description{
Returns nearest outcodes for a given longitude and latitude.
}
\examples{
\donttest{
outcode_reverse_geocoding("-3.15", "51.47")
outcode_reverse_geocoding(-3.15, 51.47)
outcode_reverse_geocoding("-3.15807731271522", "51.4799900627036")
outcode_reverse_geocoding(-3.15, 51.47, limit = 11, radius = 20000)
}

}
\seealso{
\code{\link{postcode_lookup}} for documentation.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/random_place.R
\name{random_place}
\alias{random_place}
\title{Random place}
\usage{
random_place()
}
\value{
A data frame describing a random place and all associated data.
}
\description{
Returns a random place and all associated data
}
\examples{
\donttest{
random_place()
}

}
\seealso{
\code{\link{place_lookup}} for documentation.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nearest_outcode.R
\name{nearest_outcode}
\alias{nearest_outcode}
\title{Nearest outcode}
\usage{
nearest_outcode(outcode, limit = 10, radius = 5000)
}
\arguments{
\item{outcode}{A string with a UK postcode.}

\item{limit}{An integer. Optional parameter. Limits number of postcodes matches to return.
Defaults to 10. Needs to be less than 100.}

\item{radius}{An integer. Optional parameter. Limits number of postcodes matches to return.
Defaults to 5,000m. Needs to be less than 25,000m.}
}
\value{
A list of geographical properties.
}
\description{
Returns nearest outcodes for a given outcode.
The search is based on the relative distance of the outcode centroid.
}
\examples{
\donttest{
nearest_outcode("EC1Y")
nearest_outcode("EC1Y", limit = 11)
nearest_outcode("EC1Y", limit = 11, radius = 6000)
}

}
\seealso{
\code{\link{postcode_lookup}} for documentation.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reverse_geocoding.R
\name{reverse_geocoding}
\alias{reverse_geocoding}
\title{Reverse geocoding}
\usage{
reverse_geocoding(
  longitude,
  latitude,
  limit = 10,
  radius = 100,
  wideSearch = NULL
)
}
\arguments{
\item{longitude}{A string or numeric. Needs to have at least three decimal points.}

\item{latitude}{A string or numeric. Needs to have at least three decimal points.}

\item{limit}{An integer. Limits number of postcodes matches to return. Defaults to 10. Needs to be less than 100.}

\item{radius}{An integer. Limits number of postcodes matches to return. Defaults to 100m. Needs to be less than 2,000m.}

\item{wideSearch}{TRUE or FALSE. Search up to 20km radius, but subject to a maximum of 10 results. Since lookups over a wide area can be very expensive, we've created this method to allow you choose to make the trade off between search radius and number of results. Defaults to false. When enabled, radius and limits over 10 are ignored.}
}
\value{
A list with available data.
}
\description{
Returns nearest postcodes for a given longitude and latitude.
}
\examples{
\donttest{
reverse_geocoding(0.127, 51.507)
reverse_geocoding("0.1275", "51.5073", limit = 3)
reverse_geocoding("0.1275", "51.5073", limit = 11, radius = 200)
}

}
\seealso{
\code{\link{postcode_lookup}} for documentation.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bulk_postcode_lookup.R
\name{bulk_postcode_lookup}
\alias{bulk_postcode_lookup}
\title{Bulk postcode lookup}
\usage{
bulk_postcode_lookup(postcodes)
}
\arguments{
\item{postcodes}{Accepts a list of postcodes. Accepts up to 100 postcodes.
For only one postcode use \code{\link{postcode_lookup}}.}
}
\value{
A list of length one.
}
\description{
Returns a list of matching postcodes and respective available data.
}
\examples{
\donttest{
pc_list <- list(
postcodes = c("PR3 0SG", "M45 6GN", "EX165BL")) # spaces are ignored
bulk_postcode_lookup(pc_list)
# The function needs a list of length one. This won't work:
bulk_postcode_lookup(list("PR3 0SG", "M45 6GN", "EX165BL"))
}

}
\seealso{
\code{\link{postcode_lookup}} for documentation.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/random_postcode.R
\name{random_postcode}
\alias{random_postcode}
\title{Random postcode}
\usage{
random_postcode(outcode = NULL)
}
\arguments{
\item{outcode}{A string. Filters random postcodes by outcode.
Returns null if invalid outcode. Optional.}
}
\value{
A list with a random postcode with corresponding characteristics.
}
\description{
Returns a random postcode and all available data for that postcode.
}
\examples{
\donttest{
random_postcode()
random_postcode("N1")
}

}
\seealso{
\code{\link{postcode_lookup}} for documentation.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/place_query.R
\name{place_query}
\alias{place_query}
\title{Place query}
\usage{
place_query(place, limit = 10)
}
\arguments{
\item{place}{A string. Name of a place to search for.}

\item{limit}{An integer. Limits the number of matches to return.
Defaults to 10. Needs to be less than 100.}
}
\value{
A list with available places.
}
\description{
Submit a place query and receive a complete list of places matches and associated data.
This function is similar to \code{\link{place_lookup}} but it returns a list
and allows limiting the results.
}
\examples{
\donttest{
place_query("Hills")
place_query("Hills", limit = 12)
}

}
\seealso{
\code{\link{place_lookup}} for documentation.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/postcode_validation.R
\name{postcode_validation}
\alias{postcode_validation}
\title{Postcode validation}
\usage{
postcode_validation(postcode)
}
\arguments{
\item{postcode}{A string. Valid UK postcode.}
}
\value{
A logical vector: True or False (meaning respectively valid or invalid postcode).
}
\description{
Convenience method to validate a postcode.
}
\examples{
\donttest{
postcode_validation("EC1Y 8LX") # returns TRUE
postcode_validation("XYZ") # returns FALSE
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/postcode_lookup.R
\name{postcode_lookup}
\alias{postcode_lookup}
\title{Postcode lookup}
\usage{
postcode_lookup(postcode)
}
\arguments{
\item{postcode}{A string. One valid UK postcode.
This function is case- and space-insensitive.
For more than one postcode use \code{\link{bulk_postcode_lookup}}.
For Scottish postcodes use \code{\link{scottish_postcode_lookup}}.}
}
\value{
A data frame. Returns all available data if found. Returns NAs if postcode does not exist (404).
\itemize{
\item \code{postcode} Postcode. All current ('live') postcodes within the United Kingdom,
the Channel Islands and the Isle of Man, received monthly from Royal Mail.
2, 3 or 4-character outward code, single space and 3-character inward code.
\item \code{quality} Positional Quality. Shows the status of the assigned grid reference.
\itemize{
\item 1 - within the building of the matched address closest to the postcode mean
\item 2 - as for status value 1, except by visual inspection of Landline maps (Scotland only)
\item 3 - approximate to within 50m
\item 4 - postcode unit mean (mean of matched addresses with the same postcode, but not snapped to a building)
\item 5 - imputed by ONS, by reference to surrounding postcode grid references
\item 6 - postcode sector mean, (mainly PO Boxes)
\item 8 - postcode terminated prior to Gridlink(R) initiative, last known ONS postcode grid reference1
\item 9 - no grid reference available
}
\item \code{eastings} Eastings. The Ordnance Survey postcode grid reference
Easting to 1 metre resolution; blank for postcodes in the
Channel Islands and the Isle of Man.
Grid references for postcodes in Northern Ireland relate to the Irish Grid system.
\item \code{northings} Northings. The Ordnance Survey postcode grid reference
Easting to 1 metre resolution; blank for postcodes in the
Channel Islands and the Isle of Man.
Grid references for postcodes in Northern Ireland relate to the Irish Grid system.
\item \code{country} Country. The country (i.e. one of the four constituent
countries of the United Kingdom or the Channel Islands or the Isle of Man)
to which each postcode is assigned.
\item \code{nhs_ha} Strategic Health Authority. The health area code for the postcode.
\item \code{longitude} Longitude.
The WGS84 longitude given the Postcode's national grid reference.
\item \code{latitude} Latitude.
The WGS84 latitude given the Postcode's national grid reference.
\item \code{european_electoral_region} European Electoral Region (EER).
The European Electoral Region code for each postcode.
\item \code{primary_care_trust} Primary Care Trust (PCT).
The code for the Primary Care areas in England, LHBs in Wales, CHPs in Scotland,
LCG in Northern Ireland and PHD in the Isle of Man;
there are no equivalent areas in the Channel Islands.
Care Trust/ Care Trust Plus (CT) / local health board (LHB) /
community health partnership (CHP) / local commissioning group (LCG) /
primary healthcare directorate (PHD).
\item \code{region} Region (formerly GOR). The Region code for each postcode.
The nine GORs were abolished on 1 April 2011 and are now known as 'Regions'.
They were the primary statistical subdivisions of England and
also the areas in which the Government Offices for the Regions fulfilled their role.
Each GOR covered a number of local authorities.
\item \code{lsoa} 2011 Census lower layer super output area (LSOA).
The 2011 Census lower layer SOA code for England and Wales,
SOA code for Northern Ireland and data zone code for Scotland.
\item \code{msoa} 2011 Census middle layer super output area (MSOA).
The 2011 Census middle layer SOA (MSOA) code for England and
Wales and intermediate zone for Scotland.
\item \code{incode} Incode. 3-character inward code that
is following the space in the full postcode.
\item \code{outcode} Outcode. 2, 3 or 4-character outward code.
The part of postcode before the space.
\item \code{parliamentary_constituency} Westminster Parliamentary Constituency.
The Westminster Parliamentary Constituency code for each postcode.
\item \code{admin_district} District.
The current district/unitary authority to which the postcode has been assigned.
\item \code{parish} Parish (England)/ community (Wales).
The smallest type of administrative area in England is the parish
(also known as 'civil parish'); the equivalent units in Wales are communities.
\item \code{admin_county} County. The current county to which the postcode has been assigned.
\item \code{admin_ward} Ward.
The current administrative/electoral area to which the postcode has been assigned.
\item \code{ccg} Clinical Commissioning Group. Clinical commissioning groups (CCGs)
are NHS organisations set up by the Health and Social Care Act 2012
to organise the delivery of NHS services in England.
\item \code{nuts} Nomenclature of Units for Territorial Statistics (NUTS) /
Local Administrative Units (LAU) areas.
The LAU2 code for each postcode. NUTS is a hierarchical classification of
spatial units that provides a breakdown of the European Union's
territory for producing regional statistics which are comparable across
the Union. The NUTS area classification in the United Kingdom
comprises current national administrative and electoral areas,
except in Scotland where some NUTS areas comprise whole and/or part
Local Enterprise Regions. NUTS levels 1-3 are frozen for a minimum of
three years and NUTS levels 4 and 5 are now
Local Administrative Units (LAU) levels 1 and 2 respectively.
\item \code{_code} Returns an ID or Code associated with the postcode.
Typically these are a 9 character code known as an ONS Code or GSS Code.
This is currently only available for districts, parishes, counties, CCGs, NUTS and wards.
}

See
\url{https://postcodes.io/docs} for more details.
}
\description{
Lookup a postcode.
}
\examples{
\donttest{
postcode_lookup("EC1Y8LX")
postcode_lookup("EC1Y 8LX") # spaces are ignored
postcode_lookup("DE3 5LF") # terminated postcode returns NAs
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/place_lookup.R
\name{place_lookup}
\alias{place_lookup}
\title{Place lookup}
\usage{
place_lookup(code)
}
\arguments{
\item{code}{A string. The unique identifier for places - Ordnance Survey (OSGB) code.}
}
\value{
A list with available places.
\itemize{
\item \code{code} A unique identifier that enables records to be identified easily.
The identifier will be persistent for all LocalTypes
except Section of Named Road and Section of Numbered Road.
\item \code{name_1} Name. The proper noun that applies to the real world entity.
Names that are prefixed by the definite article are not formatted
for alphabetical sorting, that is, 'The Pennines' not 'Pennines, The'.
\item \code{name_1_lang} Language of Name. The language type is only set
where more than one name exists
E.g. cym (Welsh), eng (English), gla (Scottish Gaelic).
\item \code{name_2} Name (in case of multiple languages).
The proper noun that applies to the real world entity.
Names that are prefixed by the definite article are not formatted for
alphabetical sorting, that is, 'The Pennines' not 'Pennines, The'.
\item \code{name_2_lang} Language of Name. The language type is only set where more than one name exists
E.g. cym (Welsh), eng (English), gla (Scottish Gaelic).
\item \code{local_type} The Ordnance Survey classification for the named place being represented
by the specific feature. E.g. City, Town, Village, Hamlet, Other Settlement, Suburban Area
\item \code{outcode} The postcode district, for example, SO15.
\item \code{county_unitary} Administrative Area. The name of the County (non-metropolitan or Metropolitan),
Unitary Authority or Greater London Authority administrative area that the point geometry
for feature is within or nearest to.
\item \code{county_unitary_type} Administrative Area Type. Classifies the type of administrative unit.
\item \code{district_borough} District or Borough. The name of the District, Metropolitan District or
London Borough administrative unit that the point geometry for the feature is within.
\item \code{district_borough_type} Borough Type. Classifies the type of administrative unit.
\item \code{region} The name of the European Region (was Government O ice Region)
that the point geometry for the feature is within or nearest to.
\item \code{country} The country (i.e. one of the four constituent countries of the United Kingdom
or the Channel Islands or the Isle of Man) to which each place is assigned.
\item \code{longitude} The WGS84 longitude given the Place's national grid reference.
\item \code{latitude} The WGS84 latitude given the Place's national grid reference.
\item \code{eastings} The Ordnance Survey postcode grid reference Easting to 1 metre resolution;
blank for postcodes in the Channel Islands and the Isle of Man.
\item \code{northings} The Ordnance Survey postcode grid reference Northing to 1 metre resolution;
blank for postcodes in the Channel Islands and the Isle of Man.
\item \code{min/max northings/eastings} Minimum and Maximum, Northings and Eastings.
Bounding box or Minimum Bounding Rectangle (MBR) for roads and settlements.
For Settlements and Sections of Named and Numbered Roads,
the MBR gives a representation of the extent of these features and
is not snapped to the real world extent.
}

See \url{https://postcodes.io/docs} for more details.
}
\description{
Lookup a place by code. Returns all available data if found. Returns 404 if a place does not exist.
}
\examples{
\donttest{
place_lookup("osgb4000000074544700")
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/postcode_autocomplete.R
\name{postcode_autocomplete}
\alias{postcode_autocomplete}
\title{Postcode autocomplete}
\usage{
postcode_autocomplete(postcode, limit = 10)
}
\arguments{
\item{postcode}{A string. Valid UK postcode.}

\item{limit}{An integer. Limits number of postcodes matches to return. Defaults to 10. Needs to be less than 100.}
}
\value{
A data frame with suggested postcodes.
}
\description{
Returns a data frame of matching postcodes.
}
\examples{
\donttest{
postcode_autocomplete("E1")
postcode_autocomplete("E1", limit = 11)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outward_code_lookup.R
\name{outward_code_lookup}
\alias{outward_code_lookup}
\title{Outward code lookup}
\usage{
outward_code_lookup(outcode)
}
\arguments{
\item{outcode}{A string. The outward code representing the first half of any postcode (separated by a space).}
}
\value{
The list of geographical properties.
}
\description{
Geolocation data for the centroid of the outward code specified.
}
\examples{
\donttest{
outward_code_lookup("E1")
}

}
\seealso{
\code{\link{postcode_lookup}} for documentation.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nearest_outcode_lonlat.R
\name{nearest_outcode_lonlat}
\alias{nearest_outcode_lonlat}
\title{Nearest outcodes given longitude and latitude}
\usage{
nearest_outcode_lonlat(longitude, latitude)
}
\arguments{
\item{longitude}{A string or numeric. Needs to have at least three decimal points.}

\item{latitude}{A string or numeric. Needs to have at least three decimal points.}
}
\value{
A list with available data.
}
\description{
Returns nearest outward codes for a given longitude and latitude.
The search is based on the relative distance of the outcode centroid.
}
\examples{
\donttest{
nearest_outcode_lonlat(0.127, 51.507)
}

}
\seealso{
\code{\link{postcode_lookup}} for documentation.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bulk_reverse_geocoding.R
\name{bulk_reverse_geocoding}
\alias{bulk_reverse_geocoding}
\title{Bulk reverse geocoding}
\usage{
bulk_reverse_geocoding(geolocations)
}
\arguments{
\item{geolocations}{A list containing an array of objects to geolocate.
At least two elements needed.}
}
\value{
A list with the geocoded data.
}
\description{
Returns nearest postcodes for a given longitude and latitude. Accepts up to 100 geolocations.
}
\details{
This method requires a JSON object containing an array of geolocation objects to be POSTed.
Each geolocation object accepts an optional radius (meters) and limit argument.
Both default to 100m and 10 respectively. It also accepts a wideSearch argument.
}
\examples{

\donttest{
geolocations_list <- structure(
list(
geolocations = structure(
list(
longitude = c(-3.15807731271522, -1.12935802905177),
latitude = c(51.4799900627036, 50.7186356978817),
limit = c(NA, 100L),
radius = c(NA, 500L)),
.Names = c("longitude", "latitude", "limit", "radius"),
class = "data.frame",
row.names = 1:2)),
.Names = "geolocations")

bulk_reverse_geocoding(geolocations_list)
}

}
\seealso{
\code{\link{postcode_lookup}} for documentation.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/postcode_query.R
\name{postcode_query}
\alias{postcode_query}
\title{Postcode query}
\usage{
postcode_query(postcode, limit = 10)
}
\arguments{
\item{postcode}{A string. Valid UK postcode.}

\item{limit}{An integer. Limits the number of matches to return. Defaults to 10. Needs to be less than 100.}
}
\value{
A list of geographic properties.
To return a data frame use \link[PostcodesioR]{postcode_lookup}.
}
\description{
Submit a postcode query and receive a complete list of postcode matches and all associated postcode data.
}
\examples{
\donttest{
postcode_query("EC1Y8LX")
postcode_query("EC1", limit = 11)
}

}
