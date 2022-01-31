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

# rnpn

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![cran
checks](https://cranchecks.info/badges/worst/rnpn)](https://cranchecks.info/pkgs/rnpn)
[![codecov.io](https://codecov.io/github/ropensci/rnpn/coverage.svg?branch=master)](https://codecov.io/github/ropensci/rnpn?branch=master)
[![R build
status](https://github.com/usa-npn/rnpn//workflows/R-CMD-check/badge.svg)](https://github.com/usa-npn/rnpn//actions)

`rnpn` is an R client for interacting with the USA National Phenology
Network data web services. These services include access to a rich set
of observer-contributed, point-based phenology records as well as
geospatial data products including gridded phenological model and
climatological data.

Documentation is available for the National Phenology Network [API
documentation](https://docs.google.com/document/d/1yNjupricKOAXn6tY1sI7-EwkcfwdGUZ7lxYv7fcPjO8/edit?hl=en_US),
which describes the full set of REST services this package wraps.

There is no need for an API key to grab data from the National Phenology
Network but users are required to self identify, on an honor system,
against requests that may draw upon larger datasets. Simply populate the
request\_source parameter, as necessary, with your name or the name of
your institution.

Currently there are services for writing to the database but those
endpoints do require user authentication and are not accessible through
this R wrapper. Please contact the package authors for more information
if that’s what you’re trying to do.

## Installation

This package has evolved slowly and is currently managed in a few
locations, with varying degrees of available functionality.

The original, v. 0.1, iteration of the package is available through CRAN
but has limited functionality and a number of endpoints that have been
deprecated.

CRAN version

``` r
install.packages("rnpn")
```

There’s a newer iteration of the package that includes a lot more
functionality, including the ability to access geospatial data, and is
up-to-date with the backing data services. This version of the package
is actively maintained, but is not managed through CRAN yet and as such
is more bug-prone. This version of the package must be installed through
devtools.

Development version

``` r
install.packages("devtools")
library('devtools')
devtools::install_github("usa-npn/rnpn")
```

``` r
library('rnpn')
```

## The Basics

Many of the functions to search for data require knowing the internal
unique identifiers of some of the database entities to filter the data
down efficiently. For example, if you want to search by species, then
you must know the internal identifier of the species. To get a list of
all available species use the following:

``` r
species_list <- npn_species()
```

Similarly, for phenophases:

``` r
phenophases <- npn_phenophases()
```

### Getting Observational Data

There are four main functions for accessing observational data, at
various levels of aggregation. At the most basic level you can download
the raw status and intensity data.

``` r
some_data <- npn_download_status_data(request_source='Your Name or Org Here',years=c(2015),species_id=c(35),states=c('AZ','IL'))
```

Note that through this API, data can only be filtered chronologically by
full calendar years. You can specify any number of years in each API
call. Also note that request\_source is a required parameter and should
be populated with your name or the name of the organization you
represent. All other parameters are optional but it is highly
recommended that you filter your data search further.

### Getting Geospatial Data

This package wraps around standard WCS endpoints to facilitate the
transfer of raster data. Generally, this package does not focus on
interacting with WMS services, although they are available. To get a
list of all available data layers, use the following:

``` r
layers <- npn_get_layer_details()
```

You can then use the name of the layers to select and download
geospatial data as a raster.

``` r
npn_download_geospatial(coverage_id = 'si-x:lilac_leaf_ncep_historic',date='2016-12-31',format='geotiff',output_path='./six-test-raster.tiff')
```

If you’re looking for a grid value at a specific latitude/longitude,
that is also possible.

``` r
point_value <- npn_get_point_data('si-x:lilac_leaf_ncep_historic',date='2016-12-31',lat=38.5,long=-110.7)
```

## What’s Next

Please read and review the vignettes for this package to get further
information about the full scope of functionality available.

## Meta

  - Please [report any issues or
    bugs](https://github.com/ropensci/rnpn/issues).
  - License: MIT
  - Get citation information for `rnpn` in R doing `citation(package =
    'rnpn')`
  - Please note that this package is released with a [Contributor Code
    of Conduct](https://ropensci.org/code-of-conduct/). By contributing
    to this project, you agree to abide by its terms.

[![image](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
rnpn 1.1.1 (2020-10-27)
=======================

### NEW FEATURES

* Total overhaul of the rNPN package
* Added functions for directly downloading different observation record data types
* Added additional utility and lookup type functions
* Added functions for downloading USA-NPN raster data and geospatial values by latitutde/longitude
* Deprecated the following functions: lookup_names, npn_allobssp, npn_indsatstations, npn_indspatstations, npn_species_comm, npn_species_itis, npn_species_sci, npn_stationsbystate, npn_stationswithspp



rnpn 0.1.0
==========

### NEW FEATURES

* released to CRAN
I have read and agree to the the CRAN policies at
http://cran.r-project.org/web/packages/policies.html


## Test environments

R CMD CHECK passed on local Windows 10 and Ubuntu 18 using R 4.0.2 
Also passed checks on Github Actions on macOS, Windows with R 3.6
and ubuntu 16 using R 3.5.

## R CMD check results

There were no ERRORs or WARNINGs. 

There was 1 NOTE:

Changing maintainers. This project is transitioning to a new maintainer.

## Downstream dependencies

There were no downstream dependencies.


## Other Notes

This version is a total overhaul of the previous version of the package currently available on CRAN.

Thanks! 
---
output:
  github_document
---

rnpn
========

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE,
  fig.path='inst/img/'
)
```

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/rnpn)](https://cranchecks.info/pkgs/rnpn)
[![codecov.io](https://codecov.io/github/ropensci/rnpn/coverage.svg?branch=master)](https://codecov.io/github/ropensci/rnpn?branch=master)
[![R build status](https://github.com/usa-npn/rnpn//workflows/R-CMD-check/badge.svg)](https://github.com/usa-npn/rnpn//actions)



`rnpn` is an R client for interacting with the USA National Phenology Network data web services. These services include access to a rich set of observer-contributed, point-based phenology records as well as geospatial data products including gridded phenological model and climatological data.

Documentation is available for the National Phenology Network [API documentation](https://docs.google.com/document/d/1yNjupricKOAXn6tY1sI7-EwkcfwdGUZ7lxYv7fcPjO8/edit?hl=en_US), which describes the full set of REST services this package wraps.

There is no need for an API key to grab data from the National Phenology Network but users are required to self identify, on an honor system, against requests that may draw upon larger datasets. Simply populate the request_source parameter, as necessary, with your name or the name of your institution.

Currently there are services for writing to the database but those endpoints do require user authentication and are not accessible through this R wrapper. Please contact the package authors for more information if that's what you're trying to do.

## Installation

This package has evolved slowly and is currently managed in a few locations, with varying degrees of available functionality.

The original, v. 0.1, iteration of the package is available through CRAN but has limited functionality and a number of endpoints that have been deprecated.

CRAN version

```{r eval=FALSE}
install.packages("rnpn")
```

There's a newer iteration of the package that includes a lot more functionality, including the ability to access geospatial data, and is up-to-date with the backing data services. This version of the package is actively maintained, but is not managed through CRAN yet and as such is more bug-prone. This version of the package must be installed through devtools.

Development version

```{r eval=FALSE}
install.packages("devtools")
library('devtools')
devtools::install_github("usa-npn/rnpn")
```

```{r}
library('rnpn')
```

## The Basics

Many of the functions to search for data require knowing the internal unique identifiers of some of the database entities to filter the data down efficiently. For example, if you want to search by species, then you must know the internal identifier of the species. To get a list of all available species use the following:

```{r eval=FALSE}
species_list <- npn_species()
```

Similarly, for phenophases:

```{r eval=FALSE}
phenophases <- npn_phenophases()
```

### Getting Observational Data

There are four main functions for accessing observational data, at various levels of aggregation. At the most basic level you can download the raw status and intensity data.

```{r eval=FALSE}
some_data <- npn_download_status_data(request_source='Your Name or Org Here',years=c(2015),species_id=c(35),states=c('AZ','IL'))
```

Note that through this API, data can only be filtered chronologically by full calendar years. You can specify any number of years in each API call. Also note that request_source is a required parameter and should be populated with your name or the name of the organization you represent.
All other parameters are optional but it is highly recommended that you filter your data search further.

### Getting Geospatial Data

This package wraps around standard WCS endpoints to facilitate the transfer of raster data. Generally, this package does not focus on interacting with WMS services, although they are available. To get a list of all available data layers, use the following:

```{r eval=FALSE}
layers <- npn_get_layer_details()
```
You can then use the name of the layers to select and download geospatial data as a raster.

```{r eval=FALSE}
npn_download_geospatial(coverage_id = 'si-x:lilac_leaf_ncep_historic',date='2016-12-31',format='geotiff',output_path='./six-test-raster.tiff')
```

If you're looking for a grid value at a specific latitude/longitude, that is also possible.
```{r eval=FALSE}
point_value <- npn_get_point_data('si-x:lilac_leaf_ncep_historic',date='2016-12-31',lat=38.5,long=-110.7)
```

## What's Next

Please read and review the vignettes for this package to get further information about the full scope of functionality available.

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/rnpn/issues).
* License: MIT
* Get citation information for `rnpn` in R doing `citation(package = 'rnpn')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![image](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
---
title: "VII. Putting It All Together"
author: "Lee Marsh"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteEngine{knitr::knitr}
  %\VignetteIndexEntry{07. Putting It All Together}
  %\usepackage[UTF-8]{inputenc}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(RColorBrewer)
library(rnpn)
library(rgdal)
library(raster)
```


## Combine Point and Raster Data

Observational and gridded data can be visualized or analyzed together for a variety of purposes. Users may want to identify spatial patterns in the alignment of dogwood bloom and the Spring Index bloom model. The current year's lilac leaf out observations may be compared to the 30 year average lilac sub-model of the spring index to see how well the model predicts the observations.

This example shows several data access calls to assemble observational and gridded data.



Option 1: You can add a parameter to an observational data call to additionally get a gridded layer value for each observation location/date. Note that if you don't specify which sub model of the Spring Index you want,  you will get the SI-x Average layers.

```{r eval=FALSE}
npn_download_site_phenometrics(
  request_source = 'Your Name Here', 
  years = '2013',
  num_days_quality_filter = '30', 
  species_ids = '35',
  phenophase_ids = '373', 
  download_path = 'cl_lilac_data_2013_SIxLeaf.csv',
  six_leaf_layer = TRUE,
  six_sub_model = 'lilac'
)
```

If you want to append raster data other than Spring Index, Leaf values, there's alternative boolean flags that can be set, including six_bloom_layer for Spring Index, Bloom data, and agdd_layer. Instead of TRUE or FALSE agdd_layer takes 32 or 50 and will correlate each data point with the corresponding AGDD value for the given date using either 32 or 50 base temperature.




Option 2: You can create a combined plot of observational data with modeled/raster data.

Building on the approach for accessing point data from earlier vignettes describing Individual Phenometrics and getting raster data, we can access and plot these products together. In this example, we will look at how well cloned lilac leaf out observations in 2018 are predicted by the lilac leaf sub model of the Spring Index.

### Step 1: Get the data
```{r eval=FALSE}
LilacLeaf2018<-npn_download_geospatial(
  'si-x:lilac_leaf_ncep_historic', 
  '2018-01-01', 
)


LilacLeaf2018Obs <-npn_download_individual_phenometrics(
  request_source = 'Your Name Here', 
  years = '2018',
  species_ids = '35',
  phenophase_ids = '373' 
)


```


### Step 2: Preparing the data
```{r eval=FALSE}

coords <- LilacLeaf2018Obs[ , c("longitude", "latitude")]
data <- as.data.frame(LilacLeaf2018Obs$first_yes_doy)

crs <- CRS("+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs 
                 +ellps=WGS84 +towgs84=0,0,0")

LL_spdf <- SpatialPointsDataFrame(
  coords = coords,
  data = data, 
  proj4string = crs
)

```


### Step 3: Define style options and create graph
```{r eval=FALSE}

my.palette <- brewer.pal(n=9,name="OrRd")

plot(
  LilacLeaf2018, 
  col = my.palette,
  main="2018 Observed and Predicted Lilac Leaf Out"
)

plot(
  LL_spdf,
  main="Lilac Obs",
  pch = 21,
  bg = my.palette,
  col = 'black',
  xlim=c(-125.0208,-66.47917),
  ylim=c(24.0625 ,49.9375),
  add = TRUE
)

legend(
  "bottomright", 
  legend=c("Cloned Lilac Leaf Out Observations"),
  pch = 21,
  bg = 'white',
  col = 'black',
  bty="n", 
  cex=.8
)

```
```{r, echo=FALSE, out.width = "50%", fig.pos="h"}
knitr::include_graphics("figures/7-plot.png")
---
title: "III. Individual Phenometrics"
author: "Lee Marsh"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteEngine{knitr::knitr}
  %\VignetteIndexEntry{03. Individual Phenometrics}
  %\usepackage[UTF-8]{inputenc}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(rnpn)
```

# Individual Phenometrics


While Status and Intensity data provide a direct and complete look at the observational data, some analyses rely on more synthesized output. Individual Phenometrics are derived from phenophase status data and provide estimates of phenophase onset and end dates based on the first and last "Yes" status values for organisms within a specified season of interest. Each row in this data type is comprised of values that are derived from a string of consecutive "Yes" status reports without an intervening "No" status report for a single phenophase for an individual plant or animal species at a site, called a "series". For plants, this data type provides information on the onset and end of a phenophase on an individual plant. For animals, it provides information on the onset and end of the presence of an animal species at a site. As animal presence at a site is much more likely to be interrupted by absence than the presence of a phenophase on a plant, Status and Intensity data or Site Phenometrics may be more appropriate for investigating animal phenology. However, we provide animal phenology in the same format as individual plants in this data type to allow users to readily compare individual plant phenology with animal activity.


Note that more than one series may exist for a given phenophase in an individual plant or animal species within a single growing season or year, this might occur in the case of leaf bud break followed by a killing frost and second round of breaking leaf buds. It could also occur at group sites where two or more observers are reporting on the same plant on sequential days but are not in agreement on phenophase status.


Any call for individual phenometrics requires chronological bounds, usually a calendar year, as determining onset and end depend on knowing what the time frame of interest is. If you query the services directly (without the benefit of this library) it's possible to specify arbitrary dates, in contrast this library allows you to specify a series of calendar years as input.


Here's an example of how to query the services for individual phenometrics data. Note that the overall structure and parameters are very similar to the call for status data. The biggest difference in this case is that start and end date parameters are now replaced with a 'years' array, which predictably takes a set of year values with which to query the service. 


```{r eval=FALSE}
npn_download_individual_phenometrics(
  request_source='Your Name Here', 
  years=c(2013,2014,2015,2016), 
  species_id=c(210), 
  download_path="saguaro_data_2013_2016.csv"
)
```


In this example, we're able to see individual saguaro phenology for 2013 through 2016. The results returned from the service is a tabular set of records, giving start and end date by individual saguaro plant. By default, each record contains information about the location, species, phenophase, and start and end dates.


Climate data from DayMet can also be acquired with Status & Intensity, Individual Phenometrics and Site Phenometric data types, by setting the climate_data parameter to true. In this example, we are getting colored leaves (phenophase ID is 498) data for birches, using the four birch species IDs, for 2015: 


```{r eval=FALSE}
npn_download_individual_phenometrics(
  request_source = 'Your Name Here', 
  years = c('2015'), 
  species_ids = c(97, 98, 99, 430), 
  phenophase_ids = c(498), 
  climate_data = TRUE,
  download_path = 'Betula_data_2015.csv'
)
```


To show what this looks like, we can plot the day of year of the first observation of colored leaves in birches (genus Betula) against summer Tmax.


```{r eval=FALSE}
BetulaLeaf <- read.csv(
  "Betula_data_2015.csv",
  header = TRUE, 
  na=-9999, 
  stringsAsFactors = FALSE
)

plot(
  first_yes_doy~tmax_summer, 
  data=BetulaLeaf, 
  ylab=c("Day of Year"), 
  xlab=c("Tmax Summer"),
  cex=2, 
  cex.axis=1.5, 
  cex.lab=1.5, 
  pch=21
)
```

```{r, echo=FALSE, out.width = "50%", fig.pos="h"}
knitr::include_graphics("figures/individual-phenometrics.png")
---
title: "I. Introduction to National Phenology Network Observational Data"
author: "Lee Marsh"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteEngine{knitr::knitr}
  %\VignetteIndexEntry{01. Introduction to National Phenology Network Observational Data}
  %\usepackage[UTF-8]{inputenc}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# Introduction

The USA National Phenology Network (USA-NPN) is a USGS funded organization that collects phenological observation records from volunteer and professional scientists to better understand the impact of changes in the environment on the timing of species' life cycles. The USA-NPN also provides a number of raster-based climatological data sets and phenological models. These in-situ observation and geospatial, modeled datasets are available through a number of tools and data services.

The USA-NPN R library, "rnpn", is primarily a data access service for USA-NPN data products, serving as a wrapper to the USA-NPN REST based web services. (link). This guide details how to use the library to access and work with all USA-NPN data types.


## Accessing USA-NPN Observational Data

USA-NPN Observational data are collected on the ground by citizen and professional observers following standardized protocols, using the Nature's Notebook platform. The data are available 2009 to present, and come in four formats or data types: Status & Intensity, Individual Phenometrics, Site Phenometrics and Magnitude Phenometrics. An overview of the differences is provided in the figure below, and each type is detailed in the following sections. For a complete description of the USA-NPN approach and notes for working with each data type see the [Open File Report on USA-NPN Observational Data](https://pubs.usgs.gov/of/2018/1060/ofr20181060.pdf).


In Nature's Notebook, observers register a location, and then at each location they register any number of individual plants or animal species. The expectation is that the user then takes regular observations on each individual/species at a regular interval. Phenological status is reported by yes or no answers to a series of questions, for example, "Do you see leaves?" or "Do you see active individuals?". In contrast to traditional monitoring of annual "first" events (for example, date of first leaf or first robin), this approach captures absence data when the phenophase is not occurring and repeat events. Each observation is comprised of a series of 1, 0 and -1 values, representing yes/no/uncertain for each possible phenophase for the plant on that date. To explore data in this native "Status and Intensity" format, see the vignette by the same name.

A few considerations and functions apply across all USA-NPN Observational data types. 

### Basic format for for Observational data calls

The basic format for an observational data call in the rnpn library is:


```{r eval=FALSE}
npn_download_[NAME OF DATA TYPE] (
  request_source = [NULL]
  year =  [NULL]
  species_ID = [NULL]
)                  
```

'Request source' should usually be populated with your full name or the name of the organization you represent. Species_ID is the unique identifier for all the available plants and animals in the USA-NPN database.
You can create a table of all available species and their ID numbers:

```{r eval=FALSE}
species <- npn_species()
```

Search for a species by common name from the full list:
```{r eval=FALSE}
species[species$common_name=="red maple",]
```


There are many parameters which can be set beyond these basic ones, depending on the data type, and further detailed in the other vignettes featured in this package.


### Required Parameters

Note that specifying the year(s) of interest is a required parameter.

There's also another required field, "request_source", which is a user-provided, self-identifying string. This allows the client to provide some information about who is accessing the data. Knowing who is using the data is very helpful for our staff to report the impact of the USA-NPN to the scientific community. The input provided here is entirely honor-based.


### Find stations at which a species has been observed

You can also now look up which stations have a registered plant for a particular species. In the example below, we use the species ID for red maple, which we were able to find through the npn_species() function, to find all stations with that species.

```{r eval=FALSE}
npn_stations_with_spp (3)
```

---
title: "IV. Site Phenometrics"
author: "Lee Marsh"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteEngine{knitr::knitr}
  %\VignetteIndexEntry{04. Site Phenometrics}
  %\usepackage[UTF-8]{inputenc}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(rnpn)
```

# Site Phenometrics


Site Phenometrics, derived from Individual Phenometrics, provide summary metrics of the onset and end date of phenophase activity for a species at a site. Observers are directed to create sites that represent uniform habitat and are no larger than 15 acres. For plants, this metric is calculated as an average for all individuals of a species at the site. For animals, where individuals are not tracked, this metric represents the first and last recorded appearance of the species during the season of interest. For instance, if you asked for red maple leafing data, and there was a site with three red maple trees being observed, then the data would be the average onset date for all three of those red maple trees at that site.


Here's an example of how to query the services for site phenometrics data, for cloned lilacs, breaking leaf buds, 2013. The call is very similar to the call for individual phenometrics data, however, in addition you can supply the quality control filter for the number of days between a yes record and preceding no record (also applies to the last yes and following no), for the observation to be included in the calculations. Typically this is set to 7, 14 or 30, as when downloading data using the USA-NPN Phenology Observation Portal. If you do not set this parameter, it defaults to 30 days. Note that in this example the results are stored in memory,  rather than output as a file.


```{r eval=FALSE}
LilacLeafPoints2013<-npn_download_site_phenometrics(
  request_source = 'Your Name Here', 
  years = c('2013'),
  num_days_quality_filter = '30',
  species_ids = '35',
  phenophase_ids = '373'
)

```


In this example we're able to see the date of the first observation of breaking leaf buds for cloned lilacs, averaged across individuals within sites. If any observation did not have a preceding no record within 30 days it was excluded from the calculations.


We can now plot our cloned lilac site phenometric onset data by latitude.


```{r eval=FALSE}
plot(
  mean_first_yes_doy~latitude, 
  data=LilacLeafPoints2013, 
  ylab=c("Day of Year"), 
  xlab=c("Latitude"),
  cex=2, 
  cex.axis=1.5, 
  cex.lab=1.5, 
  pch=21, 
  xlim=c(30,55), 
  ylim=c(0,200)
)
```
```{r, echo=FALSE, out.width = "50%", fig.pos="h"}
knitr::include_graphics("figures/site-phenometrics.png")
---
title: "II. Status and Intensity Data"
author: "Lee Marsh"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteEngine{knitr::knitr}
  %\VignetteIndexEntry{02. Status and Intensity Data}
  %\usepackage[UTF-8]{inputenc}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(rnpn)
```

# Status and Intensity Data

The Status and Intensity data type is the most direct presentation of the phenology data stored in the NPDb. Each row is comprised of a single record of the status (1/present/"Yes", 0/absent/"No" or -1/uncertain/"?") of a single phenophase on an individual plant or species of animal at a site on a single site visit, as well as the estimated intensity or abundance e.g., percent canopy fullness or number of individual robins observed respectively.

Retrieving this kind of data using this package is easy, and heavily parameterized. It's possible to filter data using a number of including year, geographic extent and species. In this example we get all records of bird observations in the New England states from 2018.

```{r eval=FALSE}
npn_download_status_data(
  request_source = 'Your Name Here',
  years = c('2018'),
  states = c("NY","PA","VT","MA"),
  functional_types = 'Bird'
)
```

'states' is an example of an optional parameter that allows you to filter data based on geographic location. Another example is 'functional_types' which allows you to get all available data for a group of similar species (e.g., all birds, shrubs or invasive species).
The best place to review **all** available optional filters is the autogenerated package description.

Another important optional parameter is called 'download_path'. By default requests for data from the services are returned as a data frame that gets stored in memory as a variable. In some cases, it makes more sense to save the data to file for easy and fast retrieval later. The download_path parameter allows you to specify a file path to redirect the output from the service, without having to fuss with pesky I/O operations. Additionally, requests made this way streams the data returned, so if the dataset you're working with is particularly large, it's possible to redirect the stream of data to file instead of loading it all into memory which can be useful if your environment doesn't have enough RAM to store the entire data set at once.


```{r eval=FALSE}
npn_download_status_data(
  request_source = 'Your Name Here', 
  years = c('2018'), 
  functional_types = 'Bird', 
  additional_fields = 'Site_Name', 
  download_path ='Bird_data_2018_SiteName.csv'
)
```


Using this function to get observational records is the most basic presentation of the data, and is the most robust for doing analysis, but there are a number of other products offered through the data service which provide additional value to data end users, outlined in the next vignettes.


---
title: "V. Magnitude Phenometrics"
author: "Lee Marsh"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteEngine{knitr::knitr}
  %\VignetteIndexEntry{05. Magnitude Phenometrics}
  %\usepackage[UTF-8]{inputenc}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(rnpn)
```

# Magnitude Phenometrics


Magnitude Phenometrics are a suite of eight metrics derived from Status and Intensity data. This data type provides information on the extent to which a phenophase is expressed across multiple individuals or sites, for a given set of sequential time intervals. The data user may select a weekly, bi-weekly, monthly, or custom time interval to summarize the metrics. Two metrics are available for both plants and animals, one metric is available for plants alone and five metrics are available for animals alone (table 1). Three of the five animal metrics correct animal abundance values for observer effort in time and space.


Here's an example of how to query for Magnitude Phenometrics, for the active individuals phenophase for black-capped chickadee data, in 2018. Requirements are similar to other data types. You must additionally specify the time interval by which the data should be summarized. Typically this is weekly, biweekly or monthly, as in the POP and Visualization Tool. The interval chosen in this example is 7 days.


```{r eval=FALSE}
npn_download_magnitude_phenometrics(
  request_source = 'Your Name Here', 
  years = '2018',
  period_frequency = "7",
  species_ids = '245', 
  phenophase_ids = '292', 
  download_path = 'MPM_BCC_ActInd_2018.csv'
)
```


In this example we're able to see all of the magnitude phenometric fields, including proportion_yes_records, and mean_num_animals_in-phase. See the [https://pubs.usgs.gov/of/2018/1060/ofr20181060.pdf](Open File Report on USA-NPN Observational Data) for full field descriptions.


From this dataset we can view the Proportion_Yes_Records (of all the records submitted on this species, what proportion are positive/yes records) by weekly interval:

```{r eval=FALSE}
BCC_AI<-read.csv(
  'MPM_BCC_ActInd_2018.csv', 
  header = TRUE, 
  na=-9999, 
  stringsAsFactors = FALSE
)

plot(
  BCC_AI$proportion_yes_record~as.Date(BCC_AI$start_date,"%Y-%m-%d"), 
  ylab=c("Proportion Yes Records"), 
  xlab=c("Date"),
  cex=2, 
  cex.axis=1.5, 
  cex.lab=1.5, 
  pch=21,
  xlim=as.Date(c("2018-01-01", "2018-08-01")),
  ylim=c(0,1)
)
```
```{r, echo=FALSE, out.width = "50%", fig.pos="h"}
knitr::include_graphics("figures/magnitude-phenometrics.png")


---
title: "VI. USA-NPN Geospatial Data"
author: "Lee Marsh"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteEngine{knitr::knitr}
  %\VignetteIndexEntry{06. USA-NPN Geospatial Data}
  %\usepackage[UTF-8]{inputenc}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(rnpn)
library(rgdal)
library(raster)
```

# USA-NPN Geospatial Data

USA-NPN provides phenology-relevant climate data in raster format. There are two main suites of products in this category: Accumulated Growing Degree Days and Extended Spring Indices. Accumulated Growing Degree Days and the Extended Spring Indices are both representations of accumulated temperature. As accumulated winter and spring heat drives many spring season phenological events in much of the country, these products can be used to better understand patterns in the current and historical timing of these events across the landscape. For a complete description of the USA-NPN approach and notes for working with each data type see the [Open File Report](https://pubs.usgs.gov/of/2017/1003/ofr20171003.pdf) on USA-NPN Gridded Data.

Both suites are available as:

* Current year value, with a 6-day forecast
* Current year anomaly, with a 6-day forecast
* Long-term (30 year) average
* Historical years
    + AGDD - 2016-Prior Year
    + Extended Spring Index - 1880-Prior Year
 
All of these products can be downloaded using the npn_download_geospatial call. There is a number of other products and permutations of the above listed AGDD and Spring Index products, so you can get a complete list of available layers and additional details about them including resolution, extent and the abstract/layer description. 

```{r eval=FALSE}
layers <- npn_get_layer_details()
```

The following sections describe how to parameterize calls for both AGDD and Spring Index layers. These calls result in raster data sets for the contiguous United States. 

If you are interested in how many GDDs had accumulated when the red maple in your backyard leafed out, or what day the Spring Index requirements for leaf out were met for your location, you may wish to query the layers for these values, based on location and date. There are two ways to accomplish this, using the npn_get_point_data function which works for all layers and the npn_get_AGDD_point_data function, which only works for AGDD layers and provides a more precise result.

```{r eval=FALSE}
npn_get_agdd_point_data(
 'gdd:agdd_50f', 
 '38', 
 '-90', 
 '2019-02-25'
)
```

This returns a value of 7.64098 GDD, base 50F, for the coordinates 38 north, -90 west on February 25th, 2019.


```{r eval=FALSE}
npn_get_point_data(
 'si-x:lilac_bloom_ncep', 
 '30', 
 '-90', 
 '2019-02-25'
)
```

This returns a value for lilac bloom of day 48, for the coordinates 30 north, -90 west, as of February 25th, 2019.

The above mentioned AGDD products use base temperatures of 32F or 50F and are managed through WCS services. There is also a function to get dynamic AGDD calculations based on a user defined base temperature and a number of other parameters.


```{r eval=FALSE}
custom_agdd_raster <- npn_get_custom_agdd_raster(
 method = 'double-sine',
 climate_data_source = 'NCEP',
 temp_unit = 'fahrenheit',
 start_date = '2019-01-01',
 end_date = '2019-05-10',
 base_temp = 20,
 upper_threshold = 90
)
```



## Accumulated Growing Degree Day Products

Heat accumulation is commonly used as a way of predicting the timing of phenological transitions in plants and animals, including when plants exhibit leaf out, flowering, or fruit ripening, or when insects emerge from dormancy. This is typically expressed as accumulated heat units, either Growing Degree Hours or Growing Degree Days. Growing degree day thresholds have been established for many species, and are commonly used in agriculture, horticulture, and pest management to schedule activities such as harvesting, pesticide treatment, and flower collection. The USA-NPN is currently generating Accumulated Growing Degree Days (AGDD) rasters using a January 1 start date, calculated using simple averaging. These are available calculated using two base temperatures, 32 degrees Fahrenheit (F) and 50 F. 


When querying certain layers, the underlying data is agnostic about the specific year, and in these cases it makes sense to use the day of year to request data, since that will provide a standardized result, (i.e., April 1st is day 91 in some years and day 92 in others). 

```{r eval=FALSE}
npn_download_geospatial(
 'gdd:30yr_avg_agdd_50f', 
 95
)
```

But if you're looking at a specific year, such as a current year layer, it makes sense to use a specific calendar date (formatted YYYY-MM-DD). It's also possible to save the raster directly to file instead of loading it into memory.

```{r eval=FALSE}
npn_download_geospatial(
 'gdd:agdd', 
 '2018-05-05',
 output_path='20180505-agdd-value.tiff'
)
```

In the case of the historic Spring Index layers, however, the product represents the overall outcome for the entire year, so while the year component of the date matters, the month and day do not. In this case, specify January 1 as the month and date.

```{r eval=FALSE}
npn_download_geospatial(
 "si-x:average_bloom_prism",
 "1995-01-01"
)
```

The dimension.range value, returned in the npn_get_layer_details() function, clarifies the full set of applicable dates for each layer.

Of course, it's also easy to grab raster data and load it into a visual plot as in this example, showing a map of AGDD base 50 on 2019-06-25:

```{r}
AGDDJun2019<-npn_download_geospatial(
 "gdd:agdd_50f", 
 "2019-06-25"
)
```

```{r}
plot(
 AGDDJun2019, 
 main = "AGDD base 50 on June 25th, 2019"
)
```

An important layer to know of is the 30 year average for AGDD products. This is useful for many comparative analyses. This layer takes DOY as the date input, since it's the average AGDD value for each day of year for 1981 - 2010.

```{r eval=FALSE}
average_30yr <- npn_download_geospatial(
 "gdd:30yr_avg_agdd",
 45
)
```




## Extended Spring Indices


The Extended Spring Indices are mathematical models that predict the "start of spring" (timing of first leaf or first bloom) at a particular location. These models were constructed using historical observations of the timing of first leaf and first bloom in a cloned lilac cultivar (Syringa X chinensis 'Red Rothomagensis') and two cloned honeysuckle cultivars (Lonicera tatarica L. 'Arnold Red' and Lonicera korolkowii Stapf, also known as 'Zabelii'), which were selected based on the availability of historical observations from across a wide geographic area. Primary inputs to the model are temperature and weather events, beginning January 1 of each year. The model outputs are first leaf and first bloom date for a given location. 

Data for the Spring Index is available through an enumeration of layers that represents each of the three sub-models as well as an 'average' model which represents the aggregation of the three sub-models. These layers are further enumerated by both of the represented phenophases, leaf and bloom. In the example below, first the layer representing only the Arnold Red model for 1987 is retrieved, while the second function call gets the model averaging all three of the models for the same year.

```{r eval=FALSE}
npn_download_geospatial(
 "si-x:arnoldred_bloom_prism",
 "1987-01-01"
)

average_model <- npn_download_geospatial(
 "si-x:average_bloom_prism",
 "1987-01-01"
)

```


The Spring Indices are also unique in that the algorithm has been run against the BEST climate data set, so historic data going back to 1880 is available.

```{r warning = FALSE, message=FALSE}
BESTSIxData1905 <- npn_download_geospatial(
 'si-x:average_bloom_best', 
 '1905-01-01'
)
NAvalue(BESTSIxData1905) <- -9999
```

```{r}
plot(
 BESTSIxData1905, 
 main = "Spring Index, 1905"
)

```



### Other Layers

Besides the AGDD and Spring Index layers there are a number of other useful layers available through these services, including daily temperature minimum and maximums and aggregated MODISv6 phenometrics.

The daily temperature minimum and maximum values are the underlying climate data used to generate current year AGDD and Spring Index maps. These data are generated by NOAA's National Centers for Environmental Prediction (NCEP) and are reserved through NPN's geospatial services.

```{r}
daily_max_20190505 <- npn_download_geospatial(
 'climate:tmax', 
 '2019-05-05'
)

plot(
 daily_max_20190505, 
 main = "Daily Temperature Max (C), May 5th, 2019"
)
```


The MODISv6 layers are aggregate values for remote sensing values from the MODISv6 data set, representing a subset of the following phenometrics, aggregated across 2001 - 2017: EVI Area, Mid-Greenup, Mid-Greendown. The available aggregate values for each layer are: median, TSslope, and mean absolute deviation.
This example shows the median green up value, as DOY. Note that because this layer has a fixed date, the date parameter is input as a blank string.

```{r eval=FALSE}
median_greenup <- npn_download_geospatial(
 'inca:midgup_median_nad83_02deg',
 ''
)

plot(
 median_greenup, 
 main = "MODIS Median Mid-Greenup, 2001 - 2017"
)
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npn_phenophases.R
\name{npn_abundance_categories}
\alias{npn_abundance_categories}
\title{Get Abundance Categories}
\usage{
npn_abundance_categories(...)
}
\arguments{
\item{...}{Curl options passed on to \code{\link[httr]{GET}}}
}
\description{
Gets data on all abundance/intensity categories and includes a data frame of
applicable abundance/intensity values for each category
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{npn_stationswithspp}
\alias{npn_stationswithspp}
\title{This function renamed to be consistent with other package function names.}
\usage{
npn_stationswithspp(...)
}
\description{
This function renamed to be consistent with other package function names.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npn_stations.R
\name{npn_stations}
\alias{npn_stations}
\title{Get Station Data}
\usage{
npn_stations(state_code = NULL, ...)
}
\arguments{
\item{state_code}{The postal code of the US state by which to filter
the results returned. Leave empty to get all stations.}

\item{...}{Curl options passed on to \code{\link[httr]{GET}}}
}
\value{
Stations' latitude and longitude, names, and ids.
}
\description{
Get a list of all stations, optionally filtered by state
}
\examples{
\dontrun{
npn_stations()
npn_stations('AZ')
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npn_phenophases.R
\name{npn_phenophases}
\alias{npn_phenophases}
\title{Get Phenophases}
\usage{
npn_phenophases(...)
}
\arguments{
\item{...}{Curl options passed on to \code{\link[httr]{GET}}}
}
\description{
Retrieves a complete list of all phenophases in the NPN database
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rnpn-package.R
\docType{package}
\name{rnpn-package}
\alias{rnpn-package}
\alias{rnpn}
\title{Interface to the National Phenology Network API}
\description{
This package allows for easy access to the National Phenology Network's Data API. To learn more, take a look at the vignettes.
events that occur at specific times.
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npn_phenophases.R
\name{npn_phenophases_by_species}
\alias{npn_phenophases_by_species}
\title{Get Phenophase for Species}
\usage{
npn_phenophases_by_species(species_ids, date, ...)
}
\arguments{
\item{species_ids}{List of species_ids for which to get phenophase information}

\item{date}{The applicable date for which to retrieve phenophases for the given species}

\item{...}{Curl options passed on to \code{\link[httr]{GET}}}
}
\description{
Retrieves the phenophases applicable to species for a given date. It's important to specify a date since protocols/phenophases for
any given species can change from year to year
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npn_get_species.R
\name{npn_species_types}
\alias{npn_species_types}
\title{Get Species Types}
\usage{
npn_species_types(kingdom = "Plantae", ...)
}
\arguments{
\item{kingdom}{The kingdom for which to return functional types; either 'Animalia' or 'Plantae'. Defaults to Plantae.}

\item{...}{Curl options passed on to \code{\link[httr]{GET}}}
}
\description{
Return all plant or animal functional types used in the NPN database.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npn_geoserver.R
\name{get_additional_rasters}
\alias{get_additional_rasters}
\title{Get Additional Layers}
\usage{
get_additional_rasters(data)
}
\arguments{
\item{data}{Data frame with first column named 'name' and containing the names of the layer for which to retrieve data and the second column
named 'param' and containing string representations of the time/elevation subset parameter to pass}
}
\value{
Returns a data frame containing the raster objects related to the specified layers
}
\description{
Utility function to easily take arbitrary layer name parameters as a data frame and
return the raster data from NPN Geospatial data services
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npn_data_download.R
\name{npn_download_site_phenometrics}
\alias{npn_download_site_phenometrics}
\title{Download Site Phenometrics}
\usage{
npn_download_site_phenometrics(
  request_source,
  years,
  num_days_quality_filter = "30",
  coords = NULL,
  species_ids = NULL,
  genus_ids = NULL,
  family_ids = NULL,
  order_ids = NULL,
  class_ids = NULL,
  pheno_class_ids = NULL,
  station_ids = NULL,
  species_types = NULL,
  network_ids = NULL,
  states = NULL,
  phenophase_ids = NULL,
  functional_types = NULL,
  additional_fields = NULL,
  climate_data = FALSE,
  ip_address = NULL,
  dataset_ids = NULL,
  email = NULL,
  download_path = NULL,
  six_leaf_layer = FALSE,
  six_bloom_layer = FALSE,
  agdd_layer = NULL,
  six_sub_model = NULL,
  additional_layers = NULL,
  taxonomy_aggregate = NULL,
  pheno_class_aggregate = NULL,
  wkt = NULL
)
}
\arguments{
\item{request_source}{Required field, string. Self-identify who is making requests to the data service}

\item{years}{Required field, list of strings. Specify the years to include in the search, e.g. c('2013','2014'). You must specify at least one year.}

\item{num_days_quality_filter}{Required field, defaults to 30. The integer value sets the upper limit on the number of days difference between the
first Y value and the previous N value for each individual to be included in the data aggregation.}

\item{coords}{List of float values, used to specify a bounding box as a search parameter, e.g. c ( lower_left_lat, lower_left_long,upper_right,lat,upper_right_long )}

\item{species_ids}{List of unique IDs for searching based on species, e.g. c ( 3, 34, 35 )}

\item{genus_ids}{List of unique IDs for searching based on taxonomic family, e.g. c ( 3, 34, 35 ) . This parameter will take precedence if species_ids is also set.}

\item{family_ids}{List of unique IDs for searching based on taxonomic family, e.g. c ( 3, 34, 35 ) . This parameter will take precedence if species_ids is also set.}

\item{order_ids}{List of unique IDs for searching based on taxonomic order, e.g. c ( 3, 34, 35 ) . This parameter will take precedence if species_ids or family_ids are also set.}

\item{class_ids}{List of unique IDs for searching based on taxonomic class, e.g. c ( 3, 34, 35 ) . This parameter will take precedence if species_ids, family_ids or order_ids are also set.}

\item{pheno_class_ids}{List of unique IDs for searching based on pheno class id, e.g. c (1, 5, 13)}

\item{station_ids}{List of unique IDs for searching based on site location, e.g. c ( 5, 9, ... )}

\item{species_types}{List of unique species type names for searching based on species types, e.g. c ( "Deciduous", "Evergreen" )}

\item{network_ids}{List of unique IDs for searching based on partner group/network, e.g. ( 500, 300, ... )}

\item{states}{List of US postal states to be used as search params, e.g. c ( "AZ", "IL" )}

\item{phenophase_ids}{List of unique IDs for searching based on phenophase, e.g. c ( 323, 324, ... )}

\item{functional_types}{List of unique functional type names, e.g. c ( "Birds"  )}

\item{additional_fields}{List of additional fields to be included in the search results, e.g. ( "Station_Name", "Plant_Nickname" )}

\item{climate_data}{Boolean value indicating that all climate variables should be included in additional_fields}

\item{ip_address}{Optional field, string. IP Address of user requesting data. Used for generating data reports}

\item{dataset_ids}{List of unique IDs for searching based on dataset, e.g. NEON or GRSM c(17,15)}

\item{email}{Optional field, string. Email of user requesting data.}

\item{download_path}{Optional file path to which search results should be re-directed for later use.}

\item{six_leaf_layer}{Boolean value when set to true will attempt to resolve the date of the observation to a spring index, leafing
value for the location at which the observations was taken}

\item{six_bloom_layer}{Boolean value when set to true will attempt to resolve the date of the observation to a spring index, bloom
value for the location at which the observations was taken}

\item{agdd_layer}{numeric value, accepts 32 or 50. When set, the results will attempt to resolve the date of the observation to
an AGDD value for the location; the 32 or 50 represents the base value of the AGDD value returned. All AGDD values are based on
a January 1st start date of the year in which the observation was taken.}

\item{six_sub_model}{Affects the results of the six layers returned. Can be used to specify one of three submodels used to calculate
the spring index values. Thus setting this field will change the results of six_leaf_layer and six_bloom_layer. Valid values include:
'lilac','zabelli' and 'arnoldred'. For more information see the NPN's Spring Index Maps documentation: https://www.usanpn.org/data/spring_indices}

\item{additional_layers}{Data frame with first column named 'name' and containing the names of the layer for which to retrieve data
and the second column named 'param' and containing string representations of the time/elevation subset parameter to use.
This variable can be used to append additional geospatial layer data fields to the results, such that the date of observation
in each row will resolve to a value from the specified layers, given the location of the observation.}

\item{taxonomy_aggregate}{Boolean value indicating whether to aggregate data by a taxonomic order higher than species. This will be based on the values set in family_ids, order_ids, or class_ids. If one of those three fields are not set, then this value is ignored.}

\item{pheno_class_aggregate}{Boolean value indicating whether to aggregate data by the pheno class ids as per the pheno_class_ids parameter. If the pheno_class_ids value is not set, then this parameter is ignored. This can be used in conjunction with taxonomy_aggregate and higher taxonomic level data filtering.}

\item{wkt}{WKT geometry by which filter data. Specifying a valid WKT within the contiguous US will
filter data based on the locations which fall within that WKT.}
}
\value{
Data table of all status records returned as per the search parameters. Null if output directed to file.
}
\description{
This function allows for a parameterized search of all site phenometrics records in the USA-NPN database, returning all records as per the search parameters in a
 data table. Data fetched from NPN services is returned as raw JSON before being channeled into a data table. Optionally results can be directed to an output file in
 which case raw JSON is converted to CSV and saved to file; in that case, data is also streamed to file which allows for more easily handling of the data if the search otherwise
 returns more data than can be handled at once in memory.
}
\details{
This data type includes estimates of the overall onset and end of phenophase activity for plant and animal species at a site over a user-defined time period.
 Each row provides the first and last occurrences of a given phenophase on a given species, beginning with the date of the first observed "yes" phenophase status
 record and ending with the date of the last observed "yes" record of the user-defined time period. For plant species where multiple individuals are monitored
 at the site, the date provided for "first yes" is the mean of the first "yes" records for each individual plant at the site, and the date for "last yes" is
 the mean of the last "yes" records. Note that a phenophase may have ended and restarted during the overall period of its activity at the site.
 These more fine-scale patterns can be explored in the individual phenometrics data.

 Most search parameters are optional, however, users are encouraged to supply additional search parameters to get results that are easier to work with. Request_Source
 must be provided. This is a self-identifying string, telling the service who is asking for the data or from where the request is being made. It is recommended
 you provide your name or organization name. If the call to this function is acting as an intermediary for a client, then you may also optionally provide
 a user email and/or IP address for usage data reporting later.

 Additional fields provides the ability to specify additional, non-critical fields to include in the search results. A complete list of additional fields can be found in
 the NPN service's companion documentation
 https://docs.google.com/document/d/1yNjupricKOAXn6tY1sI7-EwkcfwdGUZ7lxYv7fcPjO8/edit#heading=h.ueaexz9bczti
 Metadata on all fields can be found in the following Excel sheet:
 http://www.usanpn.org/files/metadata/site_phenometrics_datafield_descriptions.xlsx
}
\examples{
\dontrun{
#Download all saguaro data for 2013 and 2014
npn_download_site_phenometrics(
  request_source="Your Name or Org Here",
  years=c('2013','2014'),
  species_id=c(210),
  download_path="saguaro_data_2013_2014.json"
)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npn_geoserver.R
\name{resolve_six_raster}
\alias{resolve_six_raster}
\title{Resolve SIX Raster}
\usage{
resolve_six_raster(year, phenophase = "leaf", sub_model = NULL)
}
\arguments{
\item{year}{String representation of the year being requested}

\item{phenophase}{The SI-x phenophase being requested, 'leaf' or 'bloom'; defaults to 'leaf'}

\item{sub_model}{The SI-x sub model to use. Defaults to NULL (no sub-model)}
}
\value{
Returns a raster object of the appropriate SI-x layer
}
\description{
Utility function used to resolve the appropriate SI-x layer to use
based on the year being retrieved, the phenophase and sub-model being
requested.
}
\details{
If the year being requested is more than two years older than the current year
then use the prism based layers rather than the NCEP based layers.
This is because the PRISM data is not available in whole until midway through
the year after it was initially recorded. Hence, the 'safest' approach is to only
refer to the PRISM data when we knows for sure it's available in full, i.e. two years
prior.

Sub-model and phenophase on the other hand are appended to the name of the layer
to request, no special logic is present in making the decision which layer to retrieve
based on those parameters.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npn_geoserver.R
\name{npn_get_layer_details}
\alias{npn_get_layer_details}
\title{Get Geospatial Data Layer Details}
\usage{
npn_get_layer_details()
}
\value{
Data frame containing all layer details as specified in function description.
}
\description{
This function will return information about the various data layers available via the NPN's geospatial web services.
 Specifically, this function will query the NPN's GetCapabilities endpoint and parse the information on that page
 about the layers. For each layer, this function will retrieve the layer name (as to be specified elsewhere programmatically),
 the title (human readable), the abstract, which describes the data in the layer, the dimension name and dimension range for
 specifying specific date values from the layer.
}
\details{
Information about the layers can also be viewed at the getCapbilities page directly:
 https://geoserver.usanpn.org/geoserver/wms?request=GetCapabilities
}
\examples{
\dontrun{
layers <- npn_get_layer_details()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{npn_allobssp}
\alias{npn_allobssp}
\title{This function is defunct.}
\usage{
npn_allobssp(...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npn_stations.R
\name{npn_stations_by_state}
\alias{npn_stations_by_state}
\title{Get number of stations by state.}
\usage{
npn_stations_by_state(...)
}
\arguments{
\item{...}{Curl options passed on to \code{\link[httr]{GET}}}
}
\value{
Number of stations by state as a data.frame.
}
\description{
Get number of stations by state.
}
\examples{
\dontrun{
head( npn_stations_by_state() )
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npn_partner_groups.R
\name{npn_groups}
\alias{npn_groups}
\title{Get Partner Groups}
\usage{
npn_groups(use_hierarchy = FALSE, ...)
}
\arguments{
\item{use_hierarchy}{Boolean indicating whether or not the list of networks should be represented in a hierarchy. Defaults to FALSE}

\item{...}{Curl options passed on to \code{\link[httr]{GET}}}
}
\value{
List of partner groups, including ID and name
}
\description{
Returns a list of all groups participating in the NPN's data collection program. These details can be used to further filter
other service endpoints' results.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npn_dataset.R
\name{npn_datasets}
\alias{npn_datasets}
\title{Get Datasets}
\usage{
npn_datasets(...)
}
\arguments{
\item{...}{Curl options passed on to \code{\link[httr]{GET}}}
}
\value{
data.frame of datasets and their IDs
}
\description{
Returns a complete list of information about all datasets integrated into the NPN
dataset. Data can then be pulled for individual datasets using their unique IDs
}
\examples{
\dontrun{
npn_datasets()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npn_phenophases.R
\name{npn_phenophase_definitions}
\alias{npn_phenophase_definitions}
\title{Get Phenophase Definitions}
\usage{
npn_phenophase_definitions(...)
}
\arguments{
\item{...}{Curl options passed on to \code{\link[httr]{GET}}}
}
\description{
Retrieves a complete list of all phenophase definitions.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{npn_obsspbyday}
\alias{npn_obsspbyday}
\title{This function is defunct.}
\usage{
npn_obsspbyday(...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npn_data_download.R
\name{npn_download_individual_phenometrics}
\alias{npn_download_individual_phenometrics}
\title{Download Individual Phenometrics}
\usage{
npn_download_individual_phenometrics(
  request_source,
  years,
  coords = NULL,
  individual_ids = NULL,
  species_ids = NULL,
  station_ids = NULL,
  species_types = NULL,
  network_ids = NULL,
  states = NULL,
  phenophase_ids = NULL,
  functional_types = NULL,
  additional_fields = NULL,
  climate_data = FALSE,
  ip_address = NULL,
  dataset_ids = NULL,
  genus_ids = NULL,
  family_ids = NULL,
  order_ids = NULL,
  class_ids = NULL,
  pheno_class_ids = NULL,
  email = NULL,
  download_path = NULL,
  six_leaf_layer = FALSE,
  six_bloom_layer = FALSE,
  agdd_layer = NULL,
  six_sub_model = NULL,
  additional_layers = NULL,
  wkt = NULL
)
}
\arguments{
\item{request_source}{Required field, string. Self-identify who is making requests to the data service}

\item{years}{Required field, list of strings. Specify the years to include in the search, e.g. c('2013','2014'). You must specify at least one year.}

\item{coords}{List of float values, used to specify a bounding box as a search parameter, e.g. c ( lower_left_lat, lower_left_long,upper_right,lat,upper_right_long )}

\item{individual_ids}{Comma-separated string of unique IDs for individual plants/animal species by which to filter the data}

\item{species_ids}{List of unique IDs for searching based on species, e.g. c ( 3, 34, 35 )}

\item{station_ids}{List of unique IDs for searching based on site location, e.g. c ( 5, 9, ... )}

\item{species_types}{List of unique species type names for searching based on species types, e.g. c ( "Deciduous", "Evergreen" )}

\item{network_ids}{List of unique IDs for searching based on partner group/network, e.g. c( 500, 300, ... )}

\item{states}{List of US postal states to be used as search params, e.g. c ( "AZ", "IL" )}

\item{phenophase_ids}{List of unique IDs for searching based on phenophase, e.g. c ( 323, 324, ... )}

\item{functional_types}{List of unique functional type names, e.g. c ( "Birds"  )}

\item{additional_fields}{List of additional fields to be included in the search results, e.g. c ( "Station_Name", "Plant_Nickname" )}

\item{climate_data}{Boolean value indicating that all climate variables should be included in additional_fields.}

\item{ip_address}{Optional field, string. IP Address of user requesting data. Used for generating data reports}

\item{dataset_ids}{List of unique IDs for searching based on dataset, e.g. NEON or GRSM c(17,15)}

\item{genus_ids}{List of unique IDs for searching based on taxonomic family, e.g. c ( 3, 34, 35 ) . This parameter will take precedence if species_ids is also set.}

\item{family_ids}{List of unique IDs for searching based on taxonomic family, e.g. c ( 3, 34, 35 ) . This parameter will take precedence if species_ids is also set.}

\item{order_ids}{List of unique IDs for searching based on taxonomic order, e.g. c ( 3, 34, 35 ) . This parameter will take precedence if species_ids or family_ids are also set.}

\item{class_ids}{List of unique IDs for searching based on taxonomic class, e.g. c ( 3, 34, 35 ) . This parameter will take precedence if species_ids, family_ids or order_ids are also set.}

\item{pheno_class_ids}{List of unique IDs for searching based on pheno class. Note that if
both pheno_class_id and phenophase_id are provided in the same request, phenophase_id will be ignored.}

\item{email}{Optional field, string. Email of user requesting data.}

\item{download_path}{Optional file path to which search results should be re-directed for later use.}

\item{six_leaf_layer}{Boolean value when set to true will attempt to resolve the date of the observation to a spring index, leafing
value for the location at which the observations was taken}

\item{six_bloom_layer}{Boolean value when set to true will attempt to resolve the date of the observation to a spring index, bloom
value for the location at which the observations was taken}

\item{agdd_layer}{numeric value, accepts 32 or 50. When set, the results will attempt to resolve the date of the observation to
an AGDD value for the location; the 32 or 50 represents the base value of the AGDD value returned. All AGDD values are based on
a January 1st start date of the year in which the observation was taken.}

\item{six_sub_model}{Affects the results of the six layers returned. Can be used to specify one of three submodels used to calculate
the spring index values. Thus setting this field will change the results of six_leaf_layer and six_bloom_layer. Valid values include:
'lilac','zabelli' and 'arnoldred'. For more information see the NPN's Spring Index Maps documentation: https://www.usanpn.org/data/spring_indices}

\item{additional_layers}{Data frame with first column named 'name' and containing the names of the layer for which to retrieve data
and the second column named 'param' and containing string representations of the time/elevation subset parameter to use.
This variable can be used to append additional geospatial layer data fields to the results, such that the date of observation
in each row will resolve to a value from the specified layers, given the location of the observation.}

\item{wkt}{WKT geometry by which filter data. Specifying a valid WKT within the contiguous US will
filter data based on the locations which fall within that WKT.}
}
\value{
Data table of all status records returned as per the search parameters. Null if output directed to file.
}
\description{
This function allows for a parameterized search of all individual phenometrics records in the USA-NPN database, returning all records as per the search parameters in a
 data table. Data fetched from NPN services is returned as raw JSON before being channeled into a data table. Optionally results can be directed to an output file in
 which case raw JSON is converted to CSV and saved to file; in that case, data is also streamed to file which allows for more easily handling of the data if the search otherwise
 returns more data than can be handled at once in memory.
}
\details{
This data type includes estimates of the dates of phenophase onsets and ends for individual plants and for animal species at a site during a user-defined time
 period. Each row represents a series of consecutive "yes" phenophase status records, beginning with the date of the first "yes" and ending with the date of
 the last "yes", submitted for a given phenophase on a given organism. Note that more than one consecutive series for an organism may be present within a single
 growing season or year.

 Most search parameters are optional, however, users are encouraged to supply additional search parameters to get results that are easier to work with. Request_Source
 must be provided. This is a self-identifying string, telling the service who is asking for the data or from where the request is being made. It is recommended
 you provide your name or organization name. If the call to this function is acting as an intermediary for a client, then you may also optionally provide
 a user email and/or IP address for usage data reporting later.

 Additional fields provides the ability to specify additional, non-critical fields to include in the search results. A complete list of additional fields can be found in
 the NPN service's companion documentation
 https://docs.google.com/document/d/1yNjupricKOAXn6tY1sI7-EwkcfwdGUZ7lxYv7fcPjO8/edit#heading=h.7yy4i3278v7u
 Metadata on all fields can be found in the following Excel sheet:
 http://www.usanpn.org/files/metadata/individual_phenometrics_datafield_descriptions.xlsx
}
\examples{
\dontrun{
#Download all saguaro data for 2013 and 2014
npn_download_individual_phenometrics(
  request_source="Your Name or Org Here",
  years=c('2013','2014'),
  species_id=c(210),
  download_path="saguaro_data_2013_2014.json"
)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npn_geoserver.R
\name{npn_download_geospatial}
\alias{npn_download_geospatial}
\title{Download Geospatial Data}
\usage{
npn_download_geospatial(
  coverage_id,
  date,
  format = "geotiff",
  output_path = NULL
)
}
\arguments{
\item{coverage_id}{The coverage id (machine name) of the layer for which to retrieve.
Applicable values can be found via the npn_get_layer_details() function under the 'name' column.}

\item{date}{Specify the date param for the layer retrieved. This can be a calendar
date formatted YYYY-mm-dd or it could be a string integer representing day of year.
It can also be NULL in some cases. Which to use depends entirely on the layer being
requested. More information available from the npn_get_layer_details() function.}

\item{format}{The output format of the raster layer retrieved. Defaults to GeoTIFF.}

\item{output_path}{Optional value. When set, the raster will be piped to the file
path specified. When left unset, this function will return a raster object.}
}
\value{
Data frame containing all layer details as specified in function description.
}
\description{
Function for directly downloading any arbitrary Geospatial layer data from the NPN Geospatial web services.
}
\details{
Information about the layers can also be viewed at the getCapbilities page directly:
 https://geoserver.usanpn.org/geoserver/wms?request=GetCapabilities
}
\examples{
\dontrun{
ras<-npn_download_geospatial("si-x:30yr_avg_six_bloom","255")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npn_phenophases.R
\name{npn_phenophase_details}
\alias{npn_phenophase_details}
\title{Get Phenophase Details}
\usage{
npn_phenophase_details(ids, ...)
}
\arguments{
\item{ids}{List of phenophase ids for which to retrieve additional details.}

\item{...}{Curl options passed on to \code{\link[httr]{GET}}}
}
\description{
Retrieves additional details for select phenophases, including full list of applicable phenophase definition IDs and phenophase
revision notes over time
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{rnpn-defunct}
\alias{rnpn-defunct}
\title{Defunct functions in rnpn}
\description{
\itemize{
 \item \code{\link{npn_obsspbyday}}: Removed.
 \item \code{\link{npn_allobssp}}: Removed.
 \item \code{\link{npn_indspatstations}}: Removed.
 \item \code{\link{npn_indsatstations}}: Removed.
 \item \code{\link{npn_stationsbystate}}: Removed.
 \item \code{\link{npn_stationswithspp}}: Removed.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{npn_stationsbystate}
\alias{npn_stationsbystate}
\title{This function renamed to be consistent with other package function names.}
\usage{
npn_stationsbystate(...)
}
\description{
This function renamed to be consistent with other package function names.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npn_geoserver.R
\name{npn_get_custom_agdd_time_series}
\alias{npn_get_custom_agdd_time_series}
\title{Get Custom AGDD Time Series}
\usage{
npn_get_custom_agdd_time_series(
  method,
  start_date,
  end_date,
  base_temp,
  climate_data_source,
  temp_unit,
  lat,
  long,
  upper_threshold = NULL
)
}
\arguments{
\item{method}{Takes "simple" or "double-sine" as input. This is the AGDD calculation method to use for each
data point. Simple refers to simple averaging.}

\item{start_date}{Date at which to begin the AGDD calculations}

\item{end_date}{Date at which to end the AGDD calculations}

\item{base_temp}{This is the lowest temperature for each day  for it to be considered in the calculation.}

\item{climate_data_source}{Specified the climate data set to use. Takes either "PRISM" or "NCEP" as input.}

\item{temp_unit}{The unit of temperature to use in the calculation. Takes either "Fahrenheit" or "Celsius" as input.}

\item{lat}{The latitude of the location for which to calculate the time series}

\item{long}{The longitude of the location for which to calculate the time series}

\item{upper_threshold}{This parameter is only applicable for the double-sine method. This sets the highest temperature
to be considered in any given day's AGDD calculation}
}
\description{
This function takes a series of variables used in calculating AGDD and returns an AGDD time series,
based on start and end date, for a given location in the continental US.
This function leverages the USA-NPN geo web services
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npn_stations.R
\name{npn_stations_with_spp}
\alias{npn_stations_with_spp}
\title{Get Stations with Species}
\usage{
npn_stations_with_spp(speciesid, ...)
}
\arguments{
\item{speciesid}{Required. Species id numbers, from 1 to infinity, potentially,
use e.g., c(52, 53, etc.) if more than one species desired (numeric)}

\item{...}{Curl options passed on to \code{\link[httr]{GET}}}
}
\value{
Stations' latitude and longitude, names, and ids.
}
\description{
Get a list of all stations which have an individual whom is a member of a
   set of species.
}
\examples{
\dontrun{
npn_stations_with_spp(speciesid = c(52,53,54))
npn_stations_with_spp(speciesid = 53)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npn_data_download.R
\name{npn_download_status_data}
\alias{npn_download_status_data}
\title{Download Status and Intensity Records}
\usage{
npn_download_status_data(
  request_source,
  years,
  coords = NULL,
  species_ids = NULL,
  genus_ids = NULL,
  family_ids = NULL,
  order_ids = NULL,
  class_ids = NULL,
  station_ids = NULL,
  species_types = NULL,
  network_ids = NULL,
  states = NULL,
  phenophase_ids = NULL,
  functional_types = NULL,
  additional_fields = NULL,
  climate_data = FALSE,
  ip_address = NULL,
  dataset_ids = NULL,
  email = NULL,
  download_path = NULL,
  six_leaf_layer = FALSE,
  six_bloom_layer = FALSE,
  agdd_layer = NULL,
  six_sub_model = NULL,
  additional_layers = NULL,
  pheno_class_ids = NULL,
  wkt = NULL
)
}
\arguments{
\item{request_source}{Required field, string. Self-identify who is making requests to the data service}

\item{years}{Required field, list of strings. Specify the years to include in the search, e.g. c('2013','2014'). You must specify at least one year.}

\item{coords}{List of float values, used to specify a bounding box as a search parameter, e.g. c ( lower_left_lat, lower_left_long,upper_right,lat,upper_right_long )}

\item{species_ids}{List of unique IDs for searching based on species, e.g. c ( 3, 34, 35 )}

\item{genus_ids}{List of unique IDs for searching based on taxonomic family, e.g. c ( 3, 34, 35 ) . This parameter will take precedence if species_ids is also set.}

\item{family_ids}{List of unique IDs for searching based on taxonomic family, e.g. c ( 3, 34, 35 ) . This parameter will take precedence if species_ids is also set.}

\item{order_ids}{List of unique IDs for searching based on taxonomic order, e.g. c ( 3, 34, 35 ) . This parameter will take precedence if species_ids or family_ids are also set.}

\item{class_ids}{List of unique IDs for searching based on taxonomic class, e.g. c ( 3, 34, 35 ) . This parameter will take precedence if species_ids, family_ids or order_ids are also set.}

\item{station_ids}{List of unique IDs for searching based on site location, e.g. c ( 5, 9, ... )}

\item{species_types}{List of unique species type names for searching based on species types, e.g. c ( "Deciduous", "Evergreen" )}

\item{network_ids}{List of unique IDs for searching based on partner group/network, e.g. ( 500, 300, ... )}

\item{states}{List of US postal states to be used as search params, e.g. c ( "AZ", "IL" )}

\item{phenophase_ids}{List of unique IDs for searching based on phenophase, e.g. c ( 323, 324, ... )}

\item{functional_types}{List of unique functional type names, e.g. c ( "Birds"  )}

\item{additional_fields}{List of additional fields to be included in the search results, e.g. c( "Station_Name", "Plant_Nickname" )}

\item{climate_data}{Boolean value indicating that all climate variables should be included in additional_fields}

\item{ip_address}{Optional field, string. IP Address of user requesting data. Used for generating data reports}

\item{dataset_ids}{List of unique IDs for searching based on dataset, e.g. NEON or GRSM c(17,15)}

\item{email}{Optional field, string. Email of user requesting data.}

\item{download_path}{Optional file path to which search results should be re-directed for later use.}

\item{six_leaf_layer}{Boolean value when set to true will attempt to resolve the date of the observation to a spring index, leafing
value for the location at which the observations was taken}

\item{six_bloom_layer}{Boolean value when set to true will attempt to resolve the date of the observation to a spring index, bloom
value for the location at which the observations was taken}

\item{agdd_layer}{numeric value, accepts 32 or 50. When set, the results will attempt to resolve the date of the observation to
an AGDD value for the location; the 32 or 50 represents the base value of the AGDD value returned. All AGDD values are based on
a January 1st start date of the year in which the observation was taken.}

\item{six_sub_model}{Affects the results of the six layers returned. Can be used to specify one of three submodels used to calculate
the spring index values. Thus setting this field will change the results of six_leaf_layer and six_bloom_layer. Valid values include:
'lilac','zabelli' and 'arnoldred'. For more information see the NPN's Spring Index Maps documentation: https://www.usanpn.org/data/spring_indices}

\item{additional_layers}{Data frame with first column named 'name' and containing the names of the layer for which to retrieve data
and the second column named 'param' and containing string representations of the time/elevation subset parameter to use.
This variable can be used to append additional geospatial layer data fields to the results, such that the date of observation
in each row will resolve to a value from the specified layers, given the location of the observation.}

\item{pheno_class_ids}{List of unique IDs for searching based on pheno class. Note that if
both pheno_class_id and phenophase_id are provided in the same request, phenophase_id will be ignored.}

\item{wkt}{WKT geometry by which filter data. Specifying a valid WKT within the contiguous US will
filter data based on the locations which fall within that WKT.}
}
\value{
Data table of all status records returned as per the search parameters. Null if output directed to file.
}
\description{
This function allows for a parameterized search of all status records in the USA-NPN database, returning all records as per the search parameters in a data
 table. Data fetched from NPN services is returned as raw JSON before being channeled into a data table. Optionally results can be directed to an output file in
 which case the raw JSON is converted to CSV and saved to file; in that case, data is also streamed to file which allows for more easily handling of the data if the search otherwise
 returns more data than can be handled at once in memory.
}
\details{
Most search parameters are optional, however, users are encouraged to supply additional search parameters to get results that are easier to work with. Request_Source
 must be provided. This is a self-identifying string, telling the service who is asking for the data or from where the request is being made. It is recommended
 you provide your name or organization name. If the call to this function is acting as an intermediary for a client, then you may also optionally provide
 a user email and/or IP address for usage data reporting later.

 Additional fields provides the ability to specify more, non-critical fields to include in the search results. A complete list of additional fields can be found in
 the NPN service's companion documentation
 https://docs.google.com/document/d/1yNjupricKOAXn6tY1sI7-EwkcfwdGUZ7lxYv7fcPjO8/edit#heading=h.w0nctgedhaop
 Metadata on all fields can be found in the following Excel sheet:
 http://www.usanpn.org/files/metadata/status_intensity_datafield_descriptions.xlsx
}
\examples{
\dontrun{
#Download all saguaro data for 2016
npn_download_status_data(
  request_source="Your Name or Org Here",
  years=c(2016),
  species_id=c(210),
  download_path="saguaro_data_2016.json"
)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npn_data_download.R
\name{npn_get_download_url}
\alias{npn_get_download_url}
\title{Generate Download URL}
\usage{
npn_get_download_url(endpoint)
}
\arguments{
\item{endpoint}{The service point, e.g. "observations/getObservations.json?"}

\item{query_vars}{List of query params}
}
\value{
The URL, as a string
}
\description{
Utility function to create the service point URL. Base URL comes from zzz.R, endpoint is specified in the code, and query_vars should be a list of parameters.
This function will manually put those query parameters into the proper GET syntax.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npn_geoserver.R
\name{npn_get_custom_agdd_raster}
\alias{npn_get_custom_agdd_raster}
\title{Get Custom AGDD Raster Map}
\usage{
npn_get_custom_agdd_raster(
  method,
  climate_data_source,
  temp_unit,
  start_date,
  end_date,
  base_temp,
  upper_threshold = NULL
)
}
\arguments{
\item{method}{Takes "simple" or "double-sine" as input. This is the AGDD calculation method to use for each
data point. Simple refers to simple averaging.}

\item{climate_data_source}{Specified the climate data set to use. Takes either "PRISM" or "NCEP" as input.}

\item{temp_unit}{The unit of temperature to use in the calculation. Takes either "Fahrenheit" or "Celsius" as input.}

\item{start_date}{Date at which to begin the AGDD calculations}

\item{end_date}{Date at which to end the AGDD calculations}

\item{base_temp}{This is the lowest temperature for each day  for it to be considered in the calculation.}

\item{upper_threshold}{This parameter is only applicable for the double-sine method. This sets the highest temperature
to be considered in any given day's AGDD calculation}
}
\description{
This function takes a series of variables used in calculating AGDD and returns
a raster of the continental USA with each pixel representing the calculated AGDD value
based on start and end date.
This function leverages the USA-NPN geo web services
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npn_get_species.R
\name{npn_species}
\alias{npn_species}
\alias{npn_species_id}
\alias{npn_species_state}
\alias{npn_species_search}
\title{Get Species}
\usage{
npn_species(...)

npn_species_id(ids, ...)

npn_species_state(state, kingdom = NULL, ...)

npn_species_search(
  network = NULL,
  start_date = NULL,
  end_date = NULL,
  station_id = NULL,
  ...
)
}
\arguments{
\item{...}{Curl options passed on to \code{\link[httr]{GET}}}

\item{ids}{List of species ids for which to retrieve information}

\item{state}{A US postal state code to filter results}

\item{kingdom}{Filters results by taxonomic kingdom. Takes either 'Animalia' or 'Plantae'}

\item{network}{filter species based on a list of unique identifiers of NPN groups that are actually observing data on the species. Takes a list of IDs}

\item{start_date}{filter species by date observed. This sets the start date of the date range and must be used in conjunction with end_date}

\item{end_date}{filter species by date observed. This sets the end date of the date range and must be used in conjunction with start_date}

\item{station_id}{filter species by a list of unique site identifiers}
}
\value{
data.frame of species and their IDs

data.frame of the species' information
}
\description{
Returns a complete list of all species information of species represented in the NPN
database.

Returns information about a species based on the NPN's unique ID for that species

Search for species by state

Search NPN species information using a number of different parameters, which can be used in conjunction with one another, including:
 - Species on which a particular group or groups are actually collecting data
 - What species were observed in a given date range
 - What species were observed at a particular station or stations
}
\examples{
\dontrun{
npn_species()
npn_species_id(ids = 3)
}
\dontrun{
npn_species_state(state = "AZ")
npn_species_state(state = "AZ", kingdom = "Plantae")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npn_geoserver.R
\name{npn_get_agdd_point_data}
\alias{npn_get_agdd_point_data}
\title{Get AGDD Point Value}
\usage{
npn_get_agdd_point_data(layer, lat, long, date, store_data = TRUE)
}
\arguments{
\item{layer}{The name of the queried layer}

\item{lat}{The latitude of the queried point}

\item{long}{The longitude of the queried point}

\item{date}{The queried date}

\item{store_data}{Boolean value. If set TRUE then the value retrieved will be stored in a global variable named point_values for
later use}
}
\value{
Returns a numeric value of the AGDD value at the specified lat/long/date. If no value can be retrieved, then -9999 is returned.
}
\description{
This function is for requesting AGDD point values. Because the NPN has a separate
data service that can provide AGDD values which is more accurate than Geoserver
this function is ideal when requested AGDD point values.
}
\details{
As this function only works for AGDD point values, if it's necessary to retrieve point values
for other layers please try the npn_get_point_data function.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npn_geoserver.R
\name{npn_get_point_data}
\alias{npn_get_point_data}
\title{Get Point Data Value}
\usage{
npn_get_point_data(layer, lat, long, date, store_data = TRUE)
}
\arguments{
\item{layer}{The coverage id (machine name) of the layer for which to retrieve.
Applicable values can be found via the npn_get_layer_details() function under the 'name' column.}

\item{lat}{The latitude of the point}

\item{long}{The longitude of the point}

\item{date}{The date for which to get a value}

\item{store_data}{Boolean value. If set TRUE then the value retrieved will be stored in a global variable named point_values for
later use}
}
\description{
This function can get point data about any of the NPN geospatial layers.
}
\details{
Please note that this function pulls this from the NPN's WCS service so the data may not be totally precise. If
you need precise AGDD values try using the npn_get_agdd_point_data function.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{npn_indspatstations}
\alias{npn_indspatstations}
\title{This function is defunct.}
\usage{
npn_indspatstations(...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lookup_names.R
\name{npn_lookup_names}
\alias{npn_lookup_names}
\title{Species Name Lookup}
\usage{
npn_lookup_names(name, type = "genus", fuzzy = FALSE)
}
\arguments{
\item{name}{A scientific or common name}

\item{type}{One of common_name, genus, or species}

\item{fuzzy}{One of TRUE or FALSE, if FALSE, uses fuzzy search via agrep, if
FALSE, uses grep}
}
\description{
Look up species IDs by taxonomic or common name
}
\examples{
\dontrun{
npn_lookup_names(name='Pinus', type='genus')
npn_lookup_names(name='pine', type='common_name')
npn_lookup_names(name='bird', type='common_name', fuzzy=TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npn_geoserver.R
\name{npn_merge_geo_data}
\alias{npn_merge_geo_data}
\title{Merge Geo Data}
\usage{
npn_merge_geo_data(ras, col_label, df)
}
\arguments{
\item{ras}{Raster containing geospatial data}

\item{col_label}{The name of the column to append to the data frame}

\item{df}{The data frame which to append the new column of geospatial point values. For
this function to work, df must contain two columns: "longitude", and "latitude"}
}
\value{
The data frame, now appended with the new geospatial data values' column
}
\description{
Utility function to intersect point based observational data with Geospatial
data values. This will take a data frame and append a new column to it.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npn_phenophases.R
\name{npn_pheno_classes}
\alias{npn_pheno_classes}
\title{Get Pheno Classes}
\usage{
npn_pheno_classes(...)
}
\arguments{
\item{...}{Curl options passed on to \code{\link[httr]{GET}}}
}
\description{
Gets information about all pheno classes, which a higher-level order of phenophases
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{npn_indsatstations}
\alias{npn_indsatstations}
\title{This function is defunct.}
\usage{
npn_indsatstations(...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npn_phenophases.R
\name{npn_get_phenophases_for_taxon}
\alias{npn_get_phenophases_for_taxon}
\title{Get Phenophases for Taxon}
\usage{
npn_get_phenophases_for_taxon(
  family_ids = NULL,
  order_ids = NULL,
  class_ids = NULL,
  genus_ids = NULL,
  date = NULL,
  return_all = 0,
  ...
)
}
\arguments{
\item{family_ids}{List of taxonomic family ids to search for.}

\item{order_ids}{List of taxonomic order ids to search for.}

\item{class_ids}{List of taxonomic class ids to search for}

\item{genus_ids}{List of taxonomic genus ids to search for}

\item{date}{Specify the date of interest. For this function to return anything, either this value must be set of return_all must be 1.}

\item{return_all}{Takes either 0 or 1 as input and defaults to 0. For this function to return anything, either this value must be set to 1
or date must be set.}

\item{...}{Curl options passed on to \code{\link[httr]{GET}}}
}
\description{
This function gets a list of phenophases that are applicable for a provided taxonomic grouping, e.g. family, order.
Note that since a higher taxononmic order will aggregate individual species not every phenophase returned through this
function will be applicable for every species belonging to that taxonomic group.
}
\details{
It's also important to note that phenophase definitions can change for individual species over time, so there's a need
to specify either a date of interest, or to explicitly state that the function should return all phenophases that were
ever applicable for any species belonging to the specified taxonomic group.

When called, this function requires of these three parameters, exactly one of family_ids, order_ids or class_ids to be set.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npn_data_download.R
\name{npn_download_magnitude_phenometrics}
\alias{npn_download_magnitude_phenometrics}
\title{Download Magnitude Phenometrics}
\usage{
npn_download_magnitude_phenometrics(
  request_source,
  years,
  period_frequency = "30",
  coords = NULL,
  species_ids = NULL,
  genus_ids = NULL,
  family_ids = NULL,
  order_ids = NULL,
  class_ids = NULL,
  pheno_class_ids = NULL,
  station_ids = NULL,
  species_types = NULL,
  network_ids = NULL,
  states = NULL,
  phenophase_ids = NULL,
  functional_types = NULL,
  additional_fields = NULL,
  climate_data = FALSE,
  ip_address = NULL,
  dataset_ids = NULL,
  email = NULL,
  download_path = NULL,
  taxonomy_aggregate = NULL,
  pheno_class_aggregate = NULL,
  wkt = NULL
)
}
\arguments{
\item{request_source}{Required field, string. Self-identify who is making requests to the data service}

\item{years}{Required field, list of strings. Specify the years to include in the search, e.g. c('2013','2014'). You must specify at least one year.}

\item{period_frequency}{Required field, integer. The integer value specifies the number of days by which to delineate the period of time specified by the
start_date and end_date, i.e. a value of 7 will delineate the period of time weekly. Any remainder days are grouped into the final delineation.
This parameter, while typically an int, also allows for a "special" string value, "months" to be passed in. Specifying this parameter as "months" will
delineate the period of time by the calendar months regardless of how many days are in each month. Defaults to 30 if omitted.}

\item{coords}{List of float values, used to specify a bounding box as a search parameter, e.g. c ( lower_left_lat, lower_left_long,upper_right,lat,upper_right_long )}

\item{species_ids}{List of unique IDs for searching based on species, e.g. c ( 3, 34, 35 )}

\item{genus_ids}{List of unique IDs for searching based on taxonomic family, e.g. c ( 3, 34, 35 ) . This parameter will take precedence if species_ids is also set.}

\item{family_ids}{List of unique IDs for searching based on taxonomic family, e.g. c ( 3, 34, 35 ) . This parameter will take precedence if species_ids is also set.}

\item{order_ids}{List of unique IDs for searching based on taxonomic order, e.g. c ( 3, 34, 35 ) . This parameter will take precedence if species_ids or family_ids are also set.}

\item{class_ids}{List of unique IDs for searching based on taxonomic class, e.g. c ( 3, 34, 35 ) . This parameter will take precedence if species_ids, family_ids or order_ids are also set.}

\item{pheno_class_ids}{List of unique IDs for searching based on pheno class id, e.g. c (1, 5, 13)}

\item{station_ids}{List of unique IDs for searching based on site location, e.g. c ( 5, 9, ... )}

\item{species_types}{List of unique species type names for searching based on species types, e.g. c ( "Deciduous", "Evergreen" )}

\item{network_ids}{List of unique IDs for searching based on partner group/network, e.g. ( 500, 300, ... )}

\item{states}{List of US postal states to be used as search params, e.g. c ( "AZ", "IL" )}

\item{phenophase_ids}{List of unique IDs for searching based on phenophase, e.g. c ( 323, 324, ... )}

\item{functional_types}{List of unique functional type names, e.g. c ( "Birds"  )}

\item{additional_fields}{List of additional fields to be included in the search results, e.g. ( "Station_Name", "Plant_Nickname" )}

\item{climate_data}{Boolean value indicating that all climate variables should be included in additional_fields}

\item{ip_address}{Optional field, string. IP Address of user requesting data. Used for generating data reports}

\item{dataset_ids}{List of unique IDs for searching based on dataset, e.g. NEON or GRSM c(17,15)}

\item{email}{Optional field, string. Email of user requesting data.}

\item{download_path}{Optional file path to which search results should be re-directed for later use.}

\item{taxonomy_aggregate}{Boolean value indicating whether to aggregate data by a taxonomic order higher than species. This will be based on the values set in family_ids, order_ids, or class_ids. If one of those three fields are not set, then this value is ignored.}

\item{pheno_class_aggregate}{Boolean value indicating whether to aggregate data by the pheno class ids as per the pheno_class_ids parameter. If the pheno_class_ids value is not set, then this parameter is ignored. This can be used in conjunction with taxonomy_aggregate and higher taxonomic level data filtering.}

\item{wkt}{WKT geometry by which filter data. Specifying a valid WKT within the contiguous US will
filter data based on the locations which fall within that WKT.}
}
\value{
Data table of all status records returned as per the search parameters. Null if output directed to file.
}
\description{
This function allows for a parameterized search of all magnitude phenometrics in the USA-NPN database, returning all records as per the search results in a
 data table. Data fetched from NPN services is returned as raw JSON before being channeled into a data table. Optionally results can be directed to an output file in
 which case raw JSON is saved to file; in that case, data is also streamed to file which allows for more easily handling of the data if the search otherwise
 returns more data than can be handled at once in memory.
}
\details{
This data type includes various measures of the extent to which a phenophase for a plant or animal species is expressed across multiple individuals and sites
 over a user-selected set of time intervals. Each row provides up to eight calculated measures summarized weekly, bi-weekly, monthly or over a custom time interval.
 These measures include approaches to evaluate the shape of an annual activity curve, including the total number of "yes" records and the proportion of "yes"
 records relative to the total number of status records over the course of a calendar year for a region of interest. They also include several approaches for
 standardizing animal abundances by observer effort over time and space (e.g. mean active bird individuals per hour). See the Metadata window for more information.

 Most search parameters are optional, however, failing to provide even a single search parameter will return all results in the database. Request_Source
 must be provided. This is a self-identifying string, telling the service who is asking for the data or from where the request is being made. It is recommended
 you provide your name or organization name. If the call to this function is acting as an intermediary for a client, then you may also optionally provide
 a user email and/or IP address for usage data reporting later.

 Additional fields provides the ability to specify more, non-critical fields to include in the search results. A complete list of additional fields can be found in
 the NPN service's companion documentation
 https://docs.google.com/document/d/1yNjupricKOAXn6tY1sI7-EwkcfwdGUZ7lxYv7fcPjO8/edit#heading=h.df3zspopwq98
 Metadata on all fields can be found in the following Excel sheet:
 http://www.usanpn.org/files/metadata/magnitude_phenometrics_datafield_descriptions.xlsx
}
\examples{
\dontrun{
#Download book all saguaro data for 2013
npn_download_magnitude_phenometrics(
  request_source="Your Name or Org Here",
  years=c(2013),
  species_id=c(210),
  download_path="saguaro_data_2013.json"
)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npn_data_download.R
\name{npn_get_data}
\alias{npn_get_data}
\title{Download NPN Data}
\usage{
npn_get_data(
  url,
  query,
  download_path = NULL,
  always_append = FALSE,
  six_leaf_raster = NULL,
  six_bloom_raster = NULL,
  agdd_layer = NULL,
  additional_layers = NULL
)
}
\arguments{
\item{url}{The URL of the service endpoint to request data from}

\item{download_path}{String, optional file path to the file for which to output the results.}

\item{always_append}{Boolean flag. When set to true, then we always append data to the download path. This is used
in the case of npn_get_data_by_year where we're making multiple requests to the same service and aggregating all
data results in a single file. Without this flag, otherwise, each call to the service would truncate the output file.}
}
\value{
Data table of the requested data. NULL if a download_path was specified.
}
\description{
Generic utility function for querying data from the NPN data services.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npn_stations.R
\name{npn_stations_by_location}
\alias{npn_stations_by_location}
\title{Get station data based on a WKT defined geography.}
\usage{
npn_stations_by_location(wkt, ...)
}
\arguments{
\item{wkt}{Required field specifying the WKT geography to use.}

\item{...}{Curl options passed on to \code{\link[httr]{GET}}}
}
\value{
Station data as as data.frame.
}
\description{
Takes a Well-Known Text based geography as input and returns data for all
stations, including unique IDs, within that boundary.
}
\examples{
\dontrun{
head( npn_stations_by_state(wkt="POLYGON((
-110.94484396954107 32.23623109416672,-110.96166678448247 32.23594069208043,
-110.95960684795904 32.21328646993733,-110.94244071026372 32.21343170728929,
-110.93935080547857 32.23216538049456,-110.94484396954107 32.23623109416672))")
)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npn_data_download.R
\name{npn_get_common_query_vars}
\alias{npn_get_common_query_vars}
\title{Get Common Query String Variables}
\usage{
npn_get_common_query_vars(
  request_source,
  coords = NULL,
  species_ids = NULL,
  station_ids = NULL,
  species_types = NULL,
  network_ids = NULL,
  states = NULL,
  phenophase_ids = NULL,
  functional_types = NULL,
  additional_fields = NULL,
  climate_data = FALSE,
  ip_address = NULL,
  dataset_ids = NULL,
  genus_ids = NULL,
  family_ids = NULL,
  order_ids = NULL,
  class_ids = NULL,
  pheno_class_ids = NULL,
  taxonomy_aggregate = NULL,
  pheno_class_aggregate = NULL,
  wkt = NULL,
  email = NULL
)
}
\value{
List of query string variables
}
\description{
Utility function to generate a list of query string variables for requests to NPN data service points. Some parameters are basically present in all requests,
so this function helps put them together.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npn_data_download.R
\name{npn_get_data_by_year}
\alias{npn_get_data_by_year}
\title{Get Data By Year}
\usage{
npn_get_data_by_year(
  endpoint,
  query,
  years,
  download_path = NULL,
  six_leaf_layer = FALSE,
  six_bloom_layer = FALSE,
  agdd_layer = NULL,
  six_sub_model = NULL,
  additional_layers = NULL
)
}
\arguments{
\item{endpoint}{String, the endpoint to query}

\item{query}{Base query string to use. This includes all the user selected parameters but doesn't include start/end date which will be automatically generated and
added}

\item{years}{List of strings; the years for which to retrieve data. There will be one request to the service for each year}

\item{download_path}{String, optional file path to the file for which to output the results.}
}
\value{
Data table - a data table combining each requests results from the service
}
\description{
Utility function to chain multiple requests to npn_get_data for requests where data should only be retrieved on an annual basis, or otherwise automatically be
delineated in some way. Results in a data table that's a combined set of the results from each request to the data service.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.R
\name{npn_set_env}
\alias{npn_set_env}
\title{Set Environment}
\usage{
npn_set_env(env = "ops")
}
\arguments{
\item{env}{The environment to use. Should be "ops" or "dev"}
}
\description{
By default this library will call the NPN's production services
but in some cases it's preferable to access the development web services
so this function allows for manually setting the web service endpoints
to use DEV instead. Just pass in "dev" to this function to change the
endpoints to use.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npn_geoserver.R
\name{npn_check_point_cached}
\alias{npn_check_point_cached}
\title{Check Point Cached}
\usage{
npn_check_point_cached(layer, lat, long, date)
}
\arguments{
\item{layer}{The name of the queried layer}

\item{lat}{The latitude of the queried point}

\item{long}{The longitude of the queried point}

\item{date}{The queried date}
}
\value{
The value of the cell located at the specified coordinates and date if the value has been queried, otherwise NULL
}
\description{
Checks in the global variable "point values" to see if the exact data point being requested
has already been asked for and returns the value if it's already saved.
}
\keyword{internal}
