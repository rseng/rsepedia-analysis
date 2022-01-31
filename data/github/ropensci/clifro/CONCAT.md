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

Enhancing the National Climate Database with *clifro*
=====================================================

[![Build
Status](https://travis-ci.org/ropensci/clifro.svg)](https://travis-ci.org/ropensci/clifro)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/clifro)](https://cran.r-project.org/package=clifro)
[![](https://cranlogs.r-pkg.org/badges/clifro)](https://cran.r-project.org/package=clifro)
[![codecov.io](https://codecov.io/github/ropensci/clifro/coverage.svg?branch=master)](https://codecov.io/github/ropensci/clifro?branch=master)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

New Zealand’s National Climate Database,
[CliFlo](https://cliflo.niwa.co.nz/) holds data from about 6500 climate
stations, with observations dating back to 1850. CliFlo returns raw data
at ten minute, hourly, and daily frequencies. CliFlo also returns
statistical summaries, inclusive of about eighty different types of
monthly and annual statistics and six types of thirty−year normals.

The *clifro* package is designed to minimise the hassle in downloading
data from CliFlo. It does this by providing functions for the user to
log in, easily choose the appropriate datatypes and stations, and then
query the database. Once the data have been downloaded, they are stored
as specific objects in **R** with the primary aim to ensure data
visualisation and exploration is done with minimal effort and maximum
efficiency.

This package extends the functionality of
[CliFlo](https://cliflo.niwa.co.nz/) by returning stations resulting
from simultaneous searches, the ability to visualise where these climate
stations are by exporting to KML files, and elegant plotting of the
climate data. The vignettes and help files are written with the
intention that even inexperienced R users can use *clifro* easily.
Exporting the climate data from **R** is fairly easy and for more
experienced useRs, automated updating of spreadsheets or databases can
be made much easier.

Free CliFlo Subscription
------------------------

A current [CliFlo
subscription](https://cliflo.niwa.co.nz/pls/niwp/wsubform.intro) is
recommended for *clifro*, otherwise data from only one station is
available. The subscription is free and lasts for 2 years or 2,000,000
rows without renewal, which enables access to around 6,500 climate
stations around New Zealand and the Pacific.

Note this package requires internet access for connecting to the
National Climate Database web portal.

Installation in R
=================

``` r
# Install the latest CRAN release
install.packages("clifro")

# Or the latest development version
if(!require(devtools))
  install.packages("devtools")
devtools::install_github("ropensci/clifro")

# Then load the package
library(clifro)
```

Getting Started
===============

The following small example shows some of the core functionality in
*clifro*.

Where are the climate stations?
-------------------------------

We can search for climate stations anywhere in New Zealand and return
the station information in the form of a KML file. For example, we can
return all the climate stations (current and historic) in the greater
Auckland region.

``` r
all.auckland.st = cf_find_station("Auckland", search = "region", status = "all")
cf_save_kml(all.auckland.st, "all_auckland_stations")
```

![All Auckland Climate Stations](tools/README-map.png)

Note the open stations have green markers and the closed stations have
red markers.

Download and visualise public climate data
------------------------------------------

The only station available for unlimited public access to climate data
is the Reefton electronic weather station (EWS). We can download the
2014 wind and rain data and easily visualise the results very easily.

``` r
public.cfuser = cf_user()

# Choose the datatypes
daily.wind.rain.dt = cf_datatype(c(2, 3), c(1, 1), list(4, 1), c(1, NA))

# Choose the Reefton EWS station
reefton.st = cf_station()

# Send the query to CliFlo and retrieve the data
daily.datalist = cf_query(user = public.cfuser, 
                          datatype = daily.wind.rain.dt, 
                          station = reefton.st,
                          start_date = "2012-01-01 00",
                          end_date = "2013-01-01 00")
#> connecting to CliFlo...
#> reading data...
#> UserName is = public
#> Number of charged rows output = 0
#> Number of free rows output = 732
#> Total number of rows output = 732
#> Copyright NIWA 2020 Subject to NIWA's Terms and Conditions
#> See: http://clifloecd1.niwa.co.nz/pls/niwp/doc/terms.html
#> Comments to: cliflo@niwa.co.nz

# Have a look at what data is now available
daily.datalist
#> List containing clifro data frames:
#>               data      type              start                end rows
#> df 1) Surface Wind  9am only (2012-01-01  9:00) (2012-12-31  9:00)  366
#> df 2)         Rain     Daily (2012-01-01  9:00) (2012-12-31  9:00)  366

# Plot the data using default plotting methods
plot(daily.datalist)     # For the first dataframe  (Surface Wind)
```

![](tools/README-rain-wind-example-1.png)

``` r
plot(daily.datalist, 2)  # For the second dataframe (Rain)
```

![](tools/README-rain-wind-example-2.png)

For more details and reproducible examples, see the [technical
report](https://stattech.wordpress.fos.auckland.ac.nz/2015/03/25/2015-02-new-zealands-climate-data-in-r-an-introduction-to-clifro/)
for how to use *clifro*, including choosing datatypes, stations, saving
locations as KML files and easy, elegant plotting for various different
climate and weather data.

``` r
# View the clifro demo
demo(clifro)

# Read the 'Introduction to clifro' vignette
vignette("clifro")
```

Contributor Code of Conduct
===========================

The *clifro* package is released with a [contributor code of
conduct](https://github.com/ropensci/clifro/blob/master/CONDUCT.md). By
participating in this project you agree to abide by its terms.

Citation
========

``` bibtex

To cite package ‘clifro’ in publications use:

Seers B and Shears N (2015). “New Zealand's Climate Data in R - An Introduction to clifro.” The University of Auckland, Auckland, New
Zealand. <URL: https://stattech.wordpress.fos.auckland.ac.nz/2015/03/25/2015-02-new-zealands-climate-data-in-r-an-introduction-to-clifro/>.

A BibTeX entry for LaTeX users is

  @TechReport{,
    title = {New Zealand's Climate Data in R --- An Introduction to clifro},
    author = {Blake Seers and Nick Shears},
    institution = {The University of Auckland},
    address = {Auckland, New Zealand},
    year = {2015},
    url = {https://stattech.wordpress.fos.auckland.ac.nz/2015/03/25/2015-02-new-zealands-climate-data-in-r-an-introduction-to-clifro/},
  }
```

[![](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org/)
clifro 3.2.5 (18-Mar-2021)
==========
## Minor Improvements
* The `cf_query()` function now includes the `output_tz` argument that allows you to choose the output timezone as either "local" (default), "NZST", or "UTC" (issue #28).

##Bug Fixes
* Issue #27

## Dependencies
* Remove the RCurl dependency entirely from *clifro* and replace with the more modern `httr` / `rvest` / `xml2` packages.


clifro 3.2-3 (01-Sep-2020)
==========
## Bug Fixes
* Issue #21.
* Issue #22.
* Fixed broken links in documentation.
* Fixed a minor timezone bug in `cf_station`.

## Dependencies
* Remove dependency suggestion on `ggmap`. The `ggmap` library was used in the 'Working with clifro Stations' vignette to show how to plot a map of the station locations in **R**. This section of the vignette has been removed.

clifro 3.2-2 (20-Mar-2019)
==========
### Bug Fixes
* The `cf_find_station` function no longer returns an error when searching
for CliFlo stations based on proximity to a geographical coordinate (using 
the 'latlon' search) and when using a datatype (fixes issue #21).

clifro 3.2-1 (08-Jan-2019)
==========
### Bug Fixes
* Fixed issue #21

clifro 3.2-0 (25-Jul-2018)
==========
### Bug Fixes
* Fixed issue #14
* Updated links in Rd files to ensure no warnings when building package.
* `clifro` no longer tests whether or not you have Google Earth installed.

### Minor Improvements
* `clifro` has had troubles with installation on certain operating systems due 
to the `XML` package (issue #19). The `XML` and `selectr` dependencies have now 
been replaced with `xml2` and `magrittr`.
* Updated vignettes.


clifro 3.1-5 (04-Oct-2017)
==========
### Bug Fixes
* Fixed issue #14

clifro 3.1-4 (21-April-2017)
==========
### Bug Fixes
* Fixed issue #13

clifro 3.1-3 (16-March-2017)
==========
### Updates
* Reorganise package structure to include README figures to be shipped with the
package so that the upcoming version of R can compile them locally.

clifro 3.1-2 (09-January-2017)
==========
### Bug Fixes
* Fixed issue #7

clifro 3.1-0 (14-December-2016)
==========
### Minor Improvements
* Curl options can now be passed to all curl handles that are initiated by `clifro`. This means the curl options are not overwritten every time a new `clifro` function is called. Curl options are passed to `clifro` using the `cf_curl_opts` function, which is passed directly to the `RCurl::curlOptions()` function.

### Bug Fixes
* Requesting combined datatypes now works (issue #4). Note there is no default 
  plot method for this datatype as they are essentially combinations of other 
  datatypes.
* Fix bug that hung R if a datatype without any rows was requested -- Fixed issue #6

clifro 3.0-0 (10-August-2016)
==========

### Minor Improvements
* Allow expressions in legend title for windrose

### Major Bug Fixes
* HTTPS required due to a recent change in NIWA's proxy server -- Fixed Issue #3.
  As a result older versions of `clifro` don't seem to work on Windows due to an
  SSL certificate problem.

### Minor Bug Fixes
* `cf_find_station` correctly gives distances instead of longitudes

clifro 2.4-1 (15-January-2016)
==============================
### Minor Improvements
* Update citation information

clifro 2.4-0 (05-March-2015)
============================
### Bug Fixes
* Bug fixed for subsetting `cfStation` using `[`

clifro 2.2.3 (04-March-2015)
============================

* Start using NEWS to document changes to `clifro`
## Test environments
* Ubuntu Linux 16.04.6 LTS, R 4.0.0 (Travis-CI)
* Ubuntu Linux 18.04.5 LTS R 4.0.2 (local)
* Windows 10.0.19042 x64 (local)
* win-builder (devel and release)

## R CMD check results
There were no NOTEs, ERRORs, or WARNINGs.

## Downstream dependencies
I have also run R CMD CHECK on macleish, the only downstream dependency of 
clifro, without any problems.

## Resubmission

This is a resubmission. In this version I have updated the date in the 
DESCRIPTION file to match today's date.# Setup

## Platform

|setting  |value                        |
|:--------|:----------------------------|
|version  |R version 3.3.1 (2016-06-21) |
|system   |x86_64, linux-gnu            |
|ui       |RStudio (0.99.896)           |
|language |en_NZ:en                     |
|collate  |en_NZ.UTF-8                  |
|tz       |NA                           |
|date     |2016-08-10                   |

## Packages

|package      |*  |version  |date       |source         |
|:------------|:--|:--------|:----------|:--------------|
|clifro       |*  |3.0-0    |2016-08-10 |local (NA/NA)  |
|ggmap        |   |2.6.1    |2016-01-23 |CRAN (R 3.3.1) |
|ggplot2      |   |2.1.0    |2016-03-01 |CRAN (R 3.3.1) |
|knitr        |   |1.13     |2016-05-09 |CRAN (R 3.3.1) |
|lubridate    |   |1.5.6    |2016-04-06 |CRAN (R 3.3.1) |
|RColorBrewer |   |1.1-2    |2014-12-07 |CRAN (R 3.3.1) |
|RCurl        |*  |1.95-4.8 |2016-03-01 |CRAN (R 3.3.1) |
|reshape2     |   |1.4.1    |2014-12-06 |CRAN (R 3.3.1) |
|rmarkdown    |   |1.0      |2016-07-08 |CRAN (R 3.3.1) |
|scales       |   |0.4.0    |2016-02-26 |CRAN (R 3.3.1) |
|selectr      |   |0.2-3    |2014-12-24 |CRAN (R 3.3.1) |
|XML          |*  |3.98-1.4 |2016-03-01 |CRAN (R 3.3.1) |

# Check results
1 packages

## macleish (0.3.0)
Maintainer: Ben Baumer <ben.baumer@gmail.com>

0 errors | 0 warnings | 1 note 

```
checking package dependencies ... NOTE
Packages which this enhances but not available for checking:
  ‘rgdal’ ‘rgeos’
```

# Setup

## Platform

|setting  |value                        |
|:--------|:----------------------------|
|version  |R version 3.3.1 (2016-06-21) |
|system   |x86_64, linux-gnu            |
|ui       |RStudio (0.99.896)           |
|language |en_NZ:en                     |
|collate  |en_NZ.UTF-8                  |
|tz       |NA                           |
|date     |2016-08-10                   |

## Packages

|package      |*  |version  |date       |source         |
|:------------|:--|:--------|:----------|:--------------|
|clifro       |*  |3.0-0    |2016-08-10 |local (NA/NA)  |
|ggmap        |   |2.6.1    |2016-01-23 |CRAN (R 3.3.1) |
|ggplot2      |   |2.1.0    |2016-03-01 |CRAN (R 3.3.1) |
|knitr        |   |1.13     |2016-05-09 |CRAN (R 3.3.1) |
|lubridate    |   |1.5.6    |2016-04-06 |CRAN (R 3.3.1) |
|RColorBrewer |   |1.1-2    |2014-12-07 |CRAN (R 3.3.1) |
|RCurl        |*  |1.95-4.8 |2016-03-01 |CRAN (R 3.3.1) |
|reshape2     |   |1.4.1    |2014-12-06 |CRAN (R 3.3.1) |
|rmarkdown    |   |1.0      |2016-07-08 |CRAN (R 3.3.1) |
|scales       |   |0.4.0    |2016-02-26 |CRAN (R 3.3.1) |
|selectr      |   |0.2-3    |2014-12-24 |CRAN (R 3.3.1) |
|XML          |*  |3.98-1.4 |2016-03-01 |CRAN (R 3.3.1) |

# Check results
0 packages with problems


---
output:
  md_document:
    variant: markdown_github
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "tools/README-"
)
```

# Enhancing the National Climate Database with *clifro*

[![Build Status](https://travis-ci.org/ropensci/clifro.svg)](https://travis-ci.org/ropensci/clifro)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/clifro)](https://cran.r-project.org/package=clifro)
[![](https://cranlogs.r-pkg.org/badges/clifro)](https://cran.r-project.org/package=clifro)
[![codecov.io](https://codecov.io/github/ropensci/clifro/coverage.svg?branch=master)](https://codecov.io/github/ropensci/clifro?branch=master)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

New Zealand's National Climate Database, [CliFlo](https://cliflo.niwa.co.nz/) holds data from about 6500 climate stations, with observations dating back to 1850. CliFlo returns raw data at ten minute, hourly, and daily frequencies. CliFlo also returns statistical summaries, inclusive of about eighty different types of monthly and annual statistics and six types of thirty−year normals.

The *clifro* package is designed to minimise the hassle in downloading data from CliFlo. It does this by providing functions for the user to log in, easily choose the appropriate datatypes and stations, and then query the database. Once the data have been downloaded, they are stored as specific objects in **R** with the primary aim to ensure data visualisation and exploration is done with minimal effort and maximum efficiency.

This package extends the functionality of [CliFlo](https://cliflo.niwa.co.nz/) by
returning stations resulting from simultaneous searches, the ability to 
visualise where these climate stations are by exporting to KML files, and elegant
plotting of the climate data. The vignettes and help files are written with the 
intention that even inexperienced R users can use *clifro* easily. Exporting the 
climate data from **R** is fairly easy and for more experienced useRs, automated 
updating of spreadsheets or databases can be made much easier.

## Free CliFlo Subscription
A current [CliFlo subscription](https://cliflo.niwa.co.nz/pls/niwp/wsubform.intro)
is recommended for *clifro*, otherwise data from only one station is available. 
The subscription is free and lasts for 2 years or 2,000,000 rows without renewal, 
which enables access to around 6,500 climate stations around New Zealand and the 
Pacific.

Note this package requires internet access for connecting to the National 
Climate Database web portal.

# Installation in R

```{r, eval = FALSE}
# Install the latest CRAN release
install.packages("clifro")

# Or the latest development version
if(!require(devtools))
  install.packages("devtools")
devtools::install_github("ropensci/clifro")

# Then load the package
library(clifro)
```

```{r, echo = FALSE}
library(clifro)
```

# Getting Started
The following small example shows some of the core functionality in *clifro*. 

## Where are the climate stations?
We can search for climate stations anywhere in New Zealand and return the 
station information in the form of a KML file. For example, we can return all 
the climate stations (current and historic) in the greater Auckland region.

```{r, eval = FALSE}
all.auckland.st = cf_find_station("Auckland", search = "region", status = "all")
cf_save_kml(all.auckland.st, "all_auckland_stations")
```

![All Auckland Climate Stations](tools/README-map.png)

Note the open stations have green markers and the closed stations have red 
markers.

## Download and visualise public climate data

The only station available for unlimited public access to climate data is the
Reefton electronic weather station (EWS). We can download the 2014 wind and rain 
data and easily visualise the results very easily.

```{r, rain-wind-example, echo=-1}
# Create a public user
public.cfuser = cf_user()

# Choose the datatypes
daily.wind.rain.dt = cf_datatype(c(2, 3), c(1, 1), list(4, 1), c(1, NA))

# Choose the Reefton EWS station
reefton.st = cf_station()

# Send the query to CliFlo and retrieve the data
daily.datalist = cf_query(user = public.cfuser, 
                          datatype = daily.wind.rain.dt, 
                          station = reefton.st,
                          start_date = "2012-01-01 00",
                          end_date = "2013-01-01 00")

# Have a look at what data is now available
daily.datalist

# Plot the data using default plotting methods
plot(daily.datalist)     # For the first dataframe  (Surface Wind)
plot(daily.datalist, 2)  # For the second dataframe (Rain)
```

For more details and reproducible examples, see the 
[technical report](https://stattech.wordpress.fos.auckland.ac.nz/2015/03/25/2015-02-new-zealands-climate-data-in-r-an-introduction-to-clifro/)
for how to use 
*clifro*, including choosing datatypes, stations, saving locations as KML files 
and easy, elegant plotting for various different climate and weather data.

```{r, eval = FALSE}
# View the clifro demo
demo(clifro)

# Read the 'Introduction to clifro' vignette
vignette("clifro")
```

# Contributor Code of Conduct

The *clifro* package is released with a [contributor code of conduct](https://github.com/ropensci/clifro/blob/master/CONDUCT.md). By participating in this project you agree to abide by its terms.

# Citation

```bibtex

To cite package ‘clifro’ in publications use:

Seers B and Shears N (2015). “New Zealand's Climate Data in R - An Introduction to clifro.” The University of Auckland, Auckland, New
Zealand. <URL: https://stattech.wordpress.fos.auckland.ac.nz/2015/03/25/2015-02-new-zealands-climate-data-in-r-an-introduction-to-clifro/>.

A BibTeX entry for LaTeX users is

  @TechReport{,
    title = {New Zealand's Climate Data in R --- An Introduction to clifro},
    author = {Blake Seers and Nick Shears},
    institution = {The University of Auckland},
    address = {Auckland, New Zealand},
    year = {2015},
    url = {https://stattech.wordpress.fos.auckland.ac.nz/2015/03/25/2015-02-new-zealands-climate-data-in-r-an-introduction-to-clifro/},
  }
```

[![](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org/)
---
title: "Choosing a *clifro* Datatype"
author: "Blake Seers"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Choosing a *clifro* Datatype}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---
# Introduction

The `cf_datatype` function is all that is required to select `clifro` datatypes.
This function can be called without any arguments that takes the user through
interactive menus, otherwise the datatypes may be chosen programmatically if the 
menu options are known in advance. Whether the intention is to choose one 
datatype or many, this vignette details the various methods in choosing them.

# Using the menus interactively to choose a datatype

Those familiar with the cliflo [datatype selection menu](https://cliflo.niwa.co.nz/pls/niwp/wgenf.choose_datatype?cat=cat1)
will recall the myriad datatypes and options available in the National Climate 
Database. Selection of a datatype requires navigation through trees of menus,
check boxes and combo boxes. The `cf_datatype` function mimics this (tedious)
behaviour by default, i.e. when no arguments are passed to the function and the
datatypes, menus and options are all identical to (actually scraped from) 
the datatype selection menu.

## A minimal example
Let's say the datatype we are interested in is 9am surface wind in knots.

```{r, echo=FALSE}
library(clifro)
library(pander)
surfaceWind.dt = new("cfDatatype"
    , dt_name = "Wind"
    , dt_type = "Surface wind"
    , dt_sel_option_names = list("9amWind")
    , dt_sel_combo_name = "knots"
    , dt_param = structure("ls_sfw,1,2,3,4,5", .Names = "dt1")
    , dt_sel_option_params = list(structure(c("132", "knots"), .Names = c("prm4", "prm5")))
    , dt_selected_options = list(c(4, 5))
    , dt_option_length = 5
)

menu.opts = function(title, options){
  cat(paste(title, "",
              paste(seq_along(options), options, sep = ": ", 
                    collapse = "\n"), sep = "\n"))
}
```


```{r, eval=FALSE}
surfaceWind.dt = cf_datatype()

# If you prefer pointing and clicking - turn the graphics option on:
surfaceWind.dt = cf_datatype(graphics = TRUE)
```

### Daily and Hourly Observations
```{r, echo=FALSE}
menu.opts("Daily and Hourly Observations", 
          c("Combined Observations", "Wind", "Precipitation", 
                           "Temperature and Humidity", "Sunshine and Radiation", 
                           "Weather", "Pressure", "Clouds", 
                           "Evaporation / soil moisture"))
```

The first menu that appears when the above line of code is run in R is the 
'Daily and Hourly Observations'. We are interested in 'Wind', therefore we 
would type in the number of our selection (or select it using the mouse if 
`graphics = TRUE`), in this case option **2**.

### Submenu for the given datatype
```{r, echo=FALSE}
menu.opts("Wind", c("Surface wind", "Max Gust"))
```

The next menu prompts us for the type of wind we are interested in, in this case
we are interested in surface wind which is option **1**.

### Options for the given datatype
```{r, echo=FALSE}
menu.opts("Surface wind options", c("WindRun", "HlyWind", "3HlyWind", "9amWind")
          )
```

The next menu is the options for the chosen datatype, for which we may choose 
more than one. If more than one option for a given datatype is sought, options 
must be chosen one at a time. This is made possible by a menu prompting whether 
or not we would like to select another datatype every time an option is chosen.

```{r, echo=FALSE}
menu.opts("Choose another option?", c("yes", "no"))
```

We are interested only in the surface wind at 9am in this example therefore we 
don't choose another option after we choose option **4**.

### Final options
```{r, echo=FALSE}
menu.opts("Units", c("m/s", "km/hr", "knots"))
```

This final options menu is typically associated with the units of the datatype 
(although not always) and is sometimes not necessary, depending on the datatype.
For this example we do have a final menu and it prompts us for the units that 
we are interested in where we choose option **3**.

The surface wind datatype and the associated options are now saved in R as an
object called `surfaceWind.dt`.

```{r}
surfaceWind.dt
```

# Choosing a datatype without the menus
The bold numbers in the minimal example above are emphasised specifically to 
show the menu order and selections needed to choose the strength of the 9am 
surface wind in knots datatype, i.e. **2** $\rightarrow$ **1** $\rightarrow$ **4** $\rightarrow$ **3**. In 
general, if we know the selections needed for each of the four menus then we can 
choose any datatype without using the menus making datatype selection 
a lot faster and a much less tedious.

## A minimal example
To repeat our minimal example without the use of the menus we would just pass 
them as arguments to the `cf_datatype` function. These arguments are the 
selections of each of the four menus (in order) separated by a comma.

```{r, eval = FALSE}
surfaceWind.dt = cf_datatype(2, 1, 4, 3)
surfaceWind.dt
```

```{r, echo = FALSE}
surfaceWind.dt
```

## Selecting more than one option for a given datatype
Recall that we may choose more than one option at the third menu, equivalent to
the check boxes on the cliflo 
[database query form](https://cliflo.niwa.co.nz/pls/niwp/wgenf.genform1). Using 
the menu to choose more than one option is an iterative process however we can
just update our third function argument to deal with this in a more 
time-efficient manner.

```{r, echo = FALSE}
surfaceWind.dt = new("cfDatatype"
    , dt_name = "Wind"
    , dt_type = "Surface wind"
    , dt_sel_option_names = list(c("HlyWind", "9amWind"))
    , dt_sel_combo_name = "knots"
    , dt_param = structure("ls_sfw,1,2,3,4,5", .Names = "dt1")
    , dt_sel_option_params = list(structure(c("134", "132", "knots"), .Names = c("prm2", "prm4", 
"prm5")))
    , dt_selected_options = list(c(2, 4, 5))
    , dt_option_length = 5
)

rainfall.dt = new("cfDatatype"
    , dt_name = "Precipitation"
    , dt_type = "Rain (fixed periods)"
    , dt_sel_option_names = list(c("Daily ", "Hourly"))
    , dt_sel_combo_name = NA_character_
    , dt_param = structure("ls_ra,1,2,3,4", .Names = "dt1")
    , dt_sel_option_params = list(structure(c("181", "182"), .Names = c("prm1", "prm2")))
    , dt_selected_options = list(c(1, 2))
    , dt_option_length = 4
)

lightning.dt = new("cfDatatype"
    , dt_name = "Weather"
    , dt_type = "Lightning"
    , dt_sel_option_names = list("Ltng")
    , dt_sel_combo_name = NA_character_
    , dt_param = structure("ls_light,1", .Names = "dt1")
    , dt_sel_option_params = list(structure("271", .Names = "prm1"))
    , dt_selected_options = list(1)
    , dt_option_length = 1
)

temperatureExtremes.dt = new("cfDatatype"
    , dt_name = "Temperature and Humidity"
    , dt_type = "Max_min_temp"
    , dt_sel_option_names = list(c("DlyGrass", "HlyGrass"))
    , dt_sel_combo_name = NA_character_
    , dt_param = structure("ls_mxmn,1,2,3,4,5,6", .Names = "dt1")
    , dt_sel_option_params = list(structure(c("202", "204"), .Names = c("prm5", "prm6")))
    , dt_selected_options = list(c(5, 6))
    , dt_option_length = 6
)
```

```{r, eval = FALSE}
surfaceWind.dt = cf_datatype(2, 1, c(2, 4), 3)
surfaceWind.dt
```

```{r, echo = FALSE}
surfaceWind.dt
```

`surfaceWind.dt` now contains the surface wind datatype (in knots) with both 
9am wind and hourly wind. Notice how all the other function arguments remain the
same.

# Selecting multiple datatypes
Most applications involving the environmental data contained within the National
Climate Database will require selection of more than one option for more than 
one datatype. This is where the true advantages in using the `clifro` package 
become apparent.

## An extended example
Let us consider an application where we are now interested in hourly and 9am 
surface wind along with hourly and daily rainfall, hourly counts of lightning 
flashes and daily and hourly grass temperature extremes.

There are a few ways to choose all of these datatypes. Firstly, you could go 
through the menu options one by one, selecting the corresponding datatypes and 
options and saving the resulting datatypes as different R objects. A less 
laborious alternative is to create each of these datatypes without the menus.
This does of course assume we know the selections at each branch of the
[datatype selection menus](https://cliflo.niwa.co.nz/pls/niwp/wgenf.choose_datatype?cat=cat1).

```{r, eval=FALSE}
# Hourly and 9am surface wind (knots)
surfaceWind.dt = cf_datatype(2, 1, c(2, 4), 3)
surfaceWind.dt
```

```{r, echo = FALSE}
surfaceWind.dt
```

```{r, eval = FALSE}
# Hourly and daily rainfall
rainfall.dt = cf_datatype(3, 1, c(1, 2))
rainfall.dt
```

```{r, echo = FALSE}
rainfall.dt
```

```{r, eval = FALSE}
# Hourly counts of lightning flashes
lightning.dt = cf_datatype(6, 1, 1)
lightning.dt
```

```{r, echo = FALSE}
lightning.dt
```

```{r, eval = FALSE}
# Daily and hourly grass temperature extremes
temperatureExtremes.dt = cf_datatype(4, 2, c(5, 6))
temperatureExtremes.dt

# Note: only the surface wind datatype requires combo options
```

```{r, echo = FALSE}
temperatureExtremes.dt
```

This results in 4 separate objects in R containing the datatypes and their 
corresponding options. If we were wanting to submit a query using all of these 
datatypes at once, having four separate datatypes is less than optimal. The 
following table shows the options for each of the menus that we are interested
in. 

```{r, echo = FALSE, results = "asis"}
d = data.frame(Menu = c("First selection", "Second selection", 
                        "Third selection(s)", "combo box options"),
               `Surface wind` = c(2, 1, "2 & 4", 3),
               Rainfall = c(3, 1, "1 & 2", NA),
               Lightning = c(6, 1, 1, NA),
               Temperature = c(4, 2, "5 & 6", NA), check.names = FALSE)
pandoc.table(d, style = "simple")
```

We can read across the columns to see the selections that are needed to return
an R object containing the datatypes we are interested in. We can then just pass
these values into the `cf_datatype` function to return a single R object 
containing all of our datatypes and options.

```{r, echo = FALSE}
query1.dt = new("cfDatatype"
    , dt_name = c("Wind", "Precipitation", "Weather", "Temperature and Humidity"
)
    , dt_type = c("Surface wind", "Rain (fixed periods)", "Lightning", "Max_min_temp"
)
    , dt_sel_option_names = list(c("HlyWind", "9amWind"), c("Daily ", "Hourly"), "Ltng", 
    c("DlyGrass", "HlyGrass"))
    , dt_sel_combo_name = c("knots", NA, NA, NA)
    , dt_param = structure(c("ls_sfw,1,2,3,4,5", "ls_ra,6,7,8,9", "ls_light,10", 
"ls_mxmn,11,12,13,14,15,16"), .Names = c("dt1", "dt2", "dt3", 
"dt4"))
    , dt_sel_option_params = list(structure(c("134", "132", "knots"), .Names = c("prm2", "prm4", 
"prm5")), structure(c("181", "182"), .Names = c("prm6", "prm7"
)), structure("271", .Names = "prm10"), structure(c("202", "204"
), .Names = c("prm15", "prm16")))
    , dt_selected_options = list(c(2, 4, 5), c(1, 2), 1, c(5, 6))
    , dt_option_length = c(5, 4, 1, 6)
)
```

```{r, tidy = FALSE, eval = FALSE}
query1.dt = cf_datatype(c(2, 3, 6, 4), 
                        c(1, 1, 1, 2),
                        list(c(2, 4), c(1, 2), 1, c(5, 6)),
                        c(3, NA, NA, NA))
query1.dt
```

```{r, echo = FALSE}
query1.dt
```

We can also easily combine separate `cfDatatype` objects in R using the addition 
symbol `+`, to produce an identical result. This may be useful when you want 
to conduct multiple queries which include a subset of these datatypes.

```{r}
query1.dt = surfaceWind.dt + rainfall.dt + lightning.dt + 
  temperatureExtremes.dt
query1.dt
```

## Extras
```{r, eval=FALSE}
# To add another datatype using the menu:
query1.dt + cf_datatype()

# Is equivalent to:
query1.dt + cf_datatype(NA, NA, NA, NA)

# Therefore is equivalent to adding a column of NA's to the above table:
query1.dt = cf_datatype(c(2, 3, 6, 4, NA), 
                              c(1, 1, 1, 2, NA),
                              list(c(2, 4), c(1, 2), 1, c(5, 6), NA),
                              c(3, NA, NA, NA, NA))

# Half an unknown wind datatype i.e. we know first selection = 2 but nothing 
# further:
rain.dt = cf_datatype(2) # Or cf_datatype(2, NA, NA, NA)
```
---
title: "Working with *clifro* Stations"
author: "Blake Seers"
date: "`r Sys.Date()`"
output: 
  rmarkdown::html_vignette:
    fig_width: 5
    fig_height: 5
vignette: >
  %\VignetteIndexEntry{Working with *clifro* Stations}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

```{r, echo = FALSE}
library(clifro)
```
# Introduction

There are two functions available in `clifro` to create requisite `cfStation`
objects to send queries to retrieve data via `clifro`. The first one is to 
search for stations using the `cf_find_station` function as detailed in the 
[choose stations vignette][chooseStations]. The other function that creates
`cfStation` objects is the `cf_station` function that requires comma separated 
agent numbers as the only input. This vignette covers the construction of a 
`cfStation`  object via the `cf_station` function, and then shows examples of 
plotting and visualising the station's locations using KML files or within R
using the [ggmap](https://cran.r-project.org/package=ggmap)
package.

# Creating a cfStation object from agent numbers

This is the simplest method to create a `cfStation` object, simply supply the 
`cf_station` function the comma separated agent numbers. The following stations 
are (or were) located around Lake Tekapo in Canterbury, in the South Island of 
New Zealand:

1. Coal (Ski Field)
1. Macaulay (Mt Gerald)
1. South Opua
1. Mount John
1. Lake Tekapo Ews
1. Godley Peaks
1. Lilybank

```{r, eval = FALSE}
lake.tekapo.st = cf_station(12709, 35567, 39557, 4630, 24945, 4616, 4602)
lake.tekapo.st[, c("name", "agent", "start", "end", "open")]
```

```
##                      name agent      start                 end  open
## 1         Coal @ Skifield 12709 1989-02-01 2020-09-01 02:00:00  TRUE
## 2      Macaulay@Mt Gerald 35567 1990-07-04 2020-09-01 02:00:00  TRUE
## 3         Lake Tekapo Ews 24945 2003-06-18 2020-09-01 02:00:00  TRUE
## 4 South Opua @ South Opua 39557 2011-09-28 2020-09-01 02:00:00  TRUE
## 5        Lilybank Station  4602 1950-01-01 1992-09-30 00:00:00 FALSE
## 6                 Mt John  4630 1962-10-01 1988-01-01 00:00:00 FALSE
## 7    Godley Peaks, Tekapo  4616 1914-01-01 1976-06-01 00:00:00 FALSE
```

We can see that subsetting `lake.tekapo.st` acts just like a `data.frame` 
object, although it is technically a `cfStation` object. Most of the common 
`data.frame` methods work on `cfStation` objects.

## Adding more stations

To add more stations to this list the addition sign is used. Any repeated 
stations are removed and the resulting list is ordered by the end dates first 
and then by the stations' start dates.

```{r, eval = FALSE}
added.stations.st = lake.tekapo.st + 
  cf_station() + 
  cf_find_station("lighthouse", status = "all")
added.stations.st[, c("name", "agent", "start", "end", "open")]
```

```
##                       name agent      start        end  open
## 1              Reefton Ews  3925 1960-08-01 2020-09-01  TRUE
## 2          Coal @ Skifield 12709 1989-02-01 2020-09-01  TRUE
## 3       Macaulay@Mt Gerald 35567 1990-07-04 2020-09-01  TRUE
## 4          Lake Tekapo Ews 24945 2003-06-18 2020-09-01  TRUE
## 5  South Opua @ South Opua 39557 2011-09-28 2020-09-01  TRUE
## 6     Tiri Tiri Lighthouse  1401 1946-02-01 2020-08-31  TRUE
## 7  Kapoaiaia At Lighthouse 42673 1998-05-17 2020-08-31  TRUE
## 8        Orakei Lighthouse 44394 2020-05-01 2020-08-31  TRUE
## 9     Rangitoto Lighthouse 44400 2020-05-01 2020-08-31  TRUE
## 10        Lilybank Station  4602 1950-01-01 1992-09-30 FALSE
## 11                 Mt John  4630 1962-10-01 1988-01-01 FALSE
## 12   Cape Brett Lighthouse  1197 1934-11-01 1978-10-01 FALSE
## 13     Nugget Lighthouse B  5894 1975-03-01 1977-08-31 FALSE
## 14     Nugget Lighthouse A  5895 1975-03-01 1977-08-31 FALSE
## 15    Godley Peaks, Tekapo  4616 1914-01-01 1976-06-01 FALSE
## 16      Moeraki Lighthouse  5325 1935-10-01 1975-06-01 FALSE
```

The above code chunk adds the 7 stations around Lake Tekapo, the 
subscription-free reefton EWS station (`cf_station()`), and all stations presumably located 
(currently or historically) on or near a lighthouse.

Allowing multiple searches is not
currently available using the web portal, CliFlo, but the above code 
demonstrates how easy it can be in `clifro`.

# Visualising the station locations
CliFlo does not currently have any visualisation tools to aid in the selection 
of stations which can make the task of choosing geographically suitable stations
a hard one.

## Using KML files
The `cf_save_kml` functionality was introduced in the 
[choose stations vignette][chooseStations] and this function can be used on any 
`cfStation` object. To return a KML file showing all the stations within our
`added.stations.st` object we just run `cf_save_kml(added.stations.st)` in R
and the KML file is returned.

![Climate stations in the greater Auckland region.](figures/map.png)

[chooseStations]: choose-station.html
---
title: 'From CliFlo to *clifro*: An Introduction'
author: "Blake Seers"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{From CliFlo to *clifro*: An Introduction}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

```{r, echo=FALSE}
library(clifro)
```

# Introduction
The National Climate Database holds climate data from around 6,500 climate 
stations around New Zealand including some offshore and Pacific Islands. Over
600 stations are currently active and are still receiving valuable climate data.
[CliFlo](https://cliflo.niwa.co.nz/) is a web interface to the database managed 
by [NIWA](https://niwa.co.nz/), allowing users to submit queries and retrieve 
ten-minute, hourly, daily or summary data. The use of CliFlo is free given that 
the user has [subscribed](https://cliflo.niwa.co.nz/pls/niwp/wsubform.intro) and 
accepted NIWA's [terms and conditions](https://cliflo.niwa.co.nz/doc/terms.html).

The `clifro` package is designed to make CliFlo queries much simpler and
provide extensions that are currently not offered by CliFlo. The intention is 
to simplify the data extraction, manipulation, exploration and visualisation
processes and easily create publication-worthy graphics for some of the primary 
datatypes, especially for users with limited or non-existent previous R 
experience. Experienced useRs will also find this package helpful for maximising
efficiency of climate data integration with R for further analysis, modelling or 
export.

This vignette provides an introduction to the `clifro` package demonstrating the
primary functionality by way of example. For more information on any of the 
functions in the `clifro` package the user is referred to the help index for
the `clifro` package, `help(package = "clifro")`.

# Create a clifro User
As stated above, if the intention is to extract data from any station other than
Reefton Ews (subscription-free) and to maximise the potential of `clifro`, a 
valid [subscription](https://cliflo.niwa.co.nz/pls/niwp/wsubform.intro) is 
needed.

The `cf_user` function is all that is required to create a valid `clifro` user,

```{r, eval = FALSE}
me = cf_user("username", "password")
```

where `username` and `password` is substituted for the user's CliFlo 
credentials.

# Create clifro Datatypes

Once the user has been authenticated, the next step is to choose the datatypes
of interest, see the [choose datatypes vignette][chooseDatatype] for details on
choosing datatypes. For this example we are interested in daily MSL atmospheric 
pressure, minimum and maximum temperature extremes (deg C), daily rainfall (mm)
and daily surface wind.
(m/s).

```{r, eval = FALSE}
my.dts = cf_datatype(select_1 =     c(7,  4,  3,  2), 
                     select_2 =     c(1,  2,  1,  1), 
                     check_box = list(3,  1,  1,  4), 
                     combo_box =    c(NA, NA, NA, 1))
my.dts
```

```
##                      dt.name              dt.type    dt.options dt.combo
## dt1                 Pressure             Pressure      [9amMSL]         
## dt2 Temperature and Humidity         Max_min_temp [DailyMaxMin]         
## dt3            Precipitation Rain (fixed periods)      [Daily ]         
## dt4                     Wind         Surface wind     [9amWind]      m/s
```

# Create clifro Stations

The third requisite for a valid `clifro` query is the station where the
data has been collected. If the agent numbers of the required CliFlo stations
are known, the only function needed to create a clifro station `cfStation` 
object is `cf_station`. See the [choose station vignette][chooseStation]
for help with choosing stations when the agent numbers are unknown, and the 
[working with clifro stations vignette][clifrostation] for further information
and methods on `cfStation` objects.

For this example we are interested in assessing how these datatypes differ in
various parts of the country by taking a selection of stations from various
regions. These include a station from Invercargill (5814), Nelson (4241), 
Hamilton (2112) and Auckland (1962)

```{r, eval = FALSE}
my.stations = cf_station(5814, 4241, 2112, 1962)
my.stations[, 1:5]
```

```
##                name network agent      start                 end
## 1 Invercargill Aero  I68433  5814 1939-09-01 2020-08-18 02:00:00
## 2       Nelson Aero  G13222  4241 1940-07-01 2020-08-18 02:00:00
## 3     Auckland Aero  C74082  1962 1962-05-01 2020-08-18 02:00:00
## 4      Hamilton Aws  C75834  2112 1989-11-30 2020-08-18 02:00:00
```
# Retrieve the CliFlo Data

Now that we have a valid `clifro` user and the datatypes and stations of 
interest, a `clifro` query can be conducted using the `cf_query` function. We
are interested in all available data from 2012 to 2014.

```{r, eval = FALSE}
cf.datalist = cf_query(user = me, 
                       datatype = my.dts, 
                       station = my.stations, 
                       start_date = "2012-01-01 00", 
                       end_date = "2014-01-01 00")
cf.datalist
```

```
## List containing clifro data frames:
##               data      type              start                end rows
## df 1)     Pressure  9am only (2012-01-01  9:00) (2013-01-01  9:00) 1468
## df 2)      Max_min     Daily (2012-01-01  9:00) (2013-12-31  9:00) 2923
## df 3)         Rain     Daily (2012-01-01  9:00) (2013-12-31  9:00) 2923
## df 4) Surface Wind  9am only (2012-01-01  9:00) (2013-01-01  9:00) 1468
```

We can see that the pressure and surface wind data only span one year.

# Plot the CliFlo Data

There is now a list of 4 dataframes in R containing all the available data for
each of the stations and datatypes chosen above. The plotting is simply done 
with a call to `plot`, the type of plot and plotting options depends on the 
datatype. See `?'plot.cfDataList'` for details on default `clifro`
plotting. The following are examples of *some* of the plots possible with 
`clifro`, note how the optional `ggtheme` argument changes the look of the plots.

## MSL Atmospheric Pressure
This is the first dataframe in `cf.datalist`. Since the first argument passed to 
`plot` is a list of different datatypes (`cfDataList`), the second argument 
(`y`) tells the `plot` method which of the four dataframes to plot.

We could therefore simply type `plot(cf.datalist, y = 1)` and get a nice plot of
the MSL atmospheric pressure, but it is usually nice to modify the defaults 
slightly. Since the plot method returns a `ggplot` object, we can easily modify 
the plots using [ggplot2](https://ggplot2.tidyverse.org).

```{r, eval = FALSE}
# Load the ggplot2 library for element_text() and geom_smooth() functions
library(ggplot2)

# Increase the text size to 16pt and add a loess smoother with a span equal to a 
# quarter of the window
plot(cf.datalist, ggtheme = "bw", text = element_text(size = 16)) + 
  geom_smooth(method = "loess", span = 1/4)
```

![Improved MSL Atmospheric Pressure][mslAtmosPress2]

## Daily Temperature Extremes
This is the second dataframe in `cf.datalist`, therefore `y = 2`. 
These are temperature data showing the air temperature extremes at
each of the four stations, represented by a grey region in the plot. Note that 
if the average temperature were available, these would be plotted too.

```{r, eval = FALSE}
# Try a different ggtheme
plot(cf.datalist, 2, ggtheme = "linedraw")
```

![Temperature Extremes][temperature]

## Rain
This is the third dataframe in `cf.datalist`, therefore `y = 3`. Currently there 
are two possible default plots available for rainfall; with or without soil 
deficit/runoff.

```{r, eval = FALSE}
# Try yet another ggtheme
plot(cf.datalist, 3, ggtheme = "light")

# Or only plot the rainfall data
# plot(cf.datalist, 3, ggtheme = "light", include_runoff = FALSE)
```

![Rain with Soil Deficit and Runoff][rainRunoff]

## Wind

There are three types of plots available for wind data in `clifro`. The default
is to plot a windrose displaying wind speed and directions of the full time 
series at each station. The `windrose` function in `clifro` is also available 
for the user to plot their own directional data - see `?windrose`. The other two 
optional plots for wind data in `clifro` are the wind speed and wind direction 
plots. These plots display wind speed and direction patterns through time, 
adding valuable temporal information that is not portrayed in the windroses.

The wind datatype is the fourth dataframe in `cf.datalist`, therefore `y = 4`.

### Windrose

```{r, eval = FALSE}
# Defaults to windrose
plot(cf.datalist, 4, n_col = 2)
```

![Windrose][windrose]

### Wind Speeds and Directions
The other two plotting methods for wind data are the `speed_plot` and 
`direction_plot` functions to assess the temporal variability in wind (plots not 
shown).

```{r, eval = FALSE}
# Plot the wind speeds through time, choose the 'classic' ggtheme and
# allow the y-axis scales to differ for each station
speed_plot(cf.datalist, 4, ggtheme = "classic", scales = "free_y")

# Plot wind direction contours through time
direction_plot(cf.datalist, 4, n_col = 2)
```

# Data Export

```{r, eval = FALSE}
# Export the data as separate CSV files to the current working directory
for (i in seq_along(cf.datalist))
  write.csv(cf.datalist[i], 
            file = tempfile(paste0(cf.datalist[i]@dt_name, "_"), 
                            tmpdir = normalizePath("."), 
                            fileext = ".csv"),
            na = "", row.names = FALSE)

# Each dataset is saved separately here:
getwd()
```

# Summary

The primary aim of this package is to make the substantial amount of climate 
data residing within the National Climate Database more accessible and easier
to work with. The `clifro` package has many advantages over using the CliFlo 
web portal including conducting searches much more efficiently, examining the
spatial extent of the stations and enabling high quality plots to aid the data
exploration and analysis stage.

[chooseStation]: choose-station.html
[chooseDatatype]: choose-datatype.html
[clifrostation]: cfStation.html

[mslAtmosPress2]: figures/mslAtmosPress2.png "Improved MSL Atmospheric Pressure"
[temperature]: figures/temperature.png "Temperature Extremes"
[rainRunoff]: figures/rainRunoff.png "Rain with Soil Deficit and Runoff"
[windrose]: figures/windrose.png "Windrose"
---
title: "Choosing a *clifro* Station"
author: "Blake Seers"
date: "`r Sys.Date()`"
output: 
  rmarkdown::html_vignette:
    fig_width: 5
    fig_height: 5
vignette: >
  %\VignetteIndexEntry{Choosing a *clifro* Station}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

```{r, echo=FALSE}
library(clifro)
```

# Introduction

Choosing `clifro` stations is made easy with the single `cf_find_station`
function. This function is all that is required to find `clifro` stations. This
function is equivalent to conducting the same search on the 
[find stations](https://cliflo.niwa.co.nz/pls/niwp/wstn.get_stn_html) page when 
conducting a query online at CliFlo, except without some of the errors and bugs. 
This means that the searches and the types of searches possible are exactly the 
same however, `clifro` extends functionality to exploring the spatial nature of 
stations via KML files, or plotting 
directly in R. This is the main advantage in searching for stations using 
`clifro` as locating suitable stations on a map is generally the preferred 
search tool.

There are four possible types of searches:

* A search based on pattern matching the station name
* A search based on pattern matching the network ID
* A search based on region
* A search based on the vicinity of a given location

For each of these searches either all, open or closed stations may be returned
and these searches also may only return stations where given datatypes are
available. The primary goal in searching for stations is to find the 
unique station agent number required to create a `cfStation` object. This 
vignette details the various search options in `clifro` and ways to find these
requisite agent numbers, primarily by way of example.

# Ignoring datatypes
The following examples detail how to use the `cf_find_station` search function
ignoring any datatypes.

## Station name search
Both of these searches use pattern matching to find the appropriate stations. 
The station name search is useful for searching stations in certain towns or
suburbs or maybe even streets and parks. The network ID is a number that is 
assigned to the stations which makes this search useful to look up stations 
where these are known.

These searches are used when part or all of the station name or network ID is 
known. For example, consider we are looking for open stations located in Takaka, 
at the southeastern end of Golden Bay at the northern end of the South Island, 
New Zealand. The default for the `cf_find_station` function is to search *open* 
station names matching the string.

At the time of writing this, CliFlo ignores the status argument in the name and 
network ID search whereas `clifro` does not. Searching open stations with the
station name matching "takaka" on CliFlo will return these stations.

```{r, eval = FALSE}
# Equivalent to searching for status = "open" on CliFro
# Note the search string is not case sensitive
cf_find_station("takaka", status = "all")
```

```{r, echo = FALSE}
takaka.df = structure(list(name = c("Takaka, Kotinga Road", "Riwaka At Takaka Hill", 
"Takaka Pohara", "Takaka At Harwoods", "Takaka At Kotinga", "Takaka @ Canaan", 
"Upper Takaka 2", "Takaka Ews", "Takaka Aero Raws", "Takaka, Kotinga 2", 
"Upper Takaka", "Takaka,Patons Rock", "Takaka,Kotinga 1", "Takaka Aero", 
"Takaka Hill", "Takaka,Bu Bu", "Takaka"), network = c("F02882", 
"O12090", "F02884", "F15292", "F15291", "F0299A", "F12083", "F02885", 
"O00957", "F02883", "F12082", "F02772", "F02971", "F02871", "F12081", 
"F02872", "F02881"), agent = c(3788L, 44046L, 3790L, 44050L, 
44051L, 44072L, 11519L, 23849L, 41196L, 3789L, 7316L, 3779L, 
3794L, 3785L, 3833L, 3786L, 3787L), start = structure(c(18273600, 
316263600, 520516800, 570020400, 704030400, 760014000, 805464000, 
1020081600, 1439294400, 502110000, 688820400, -7992000, -255182400, 
-1046692800, -704894400, -1159875000, -2082886200), class = c("POSIXct", 
"POSIXt"), tzone = "NZ"), end = structure(c(1597665600, 1597665600, 
1597665600, 1597665600, 1597665600, 1597665600, 1597665600, 1597665600, 
1597665600, 1341057600, 720442800, 157719600, 49809600, 7732800, 
-320932800, -760190400, -1333452600), class = c("POSIXct", "POSIXt"
), tzone = "NZ"), open = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, 
TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
FALSE), distance = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
NA, NA, NA, NA, NA, NA, NA), lat = c(-40.872, -41.03192, -40.845, 
-41.03094, -40.87068, -40.93987, -41.01516, -40.86364, -40.81531, 
-40.882, -41.051, -40.789, -40.9, -40.816, -41.017, -40.85, -40.817
), lon = c(172.809, 172.84439, 172.867, 172.79802, 172.808, 172.90821, 
172.82582, 172.80568, 172.7765, 172.801, 172.833, 172.757, 172.775, 
172.772, 172.867, 172.733, 172.8)), class = "data.frame", row.names = c(NA, 
-17L))

new("cfStation", takaka.df)
```

This shows that 8 of these 17 stations are closed. The search in `clifro` does 
not ignore the station status.

```{r, eval = FALSE}
cf_find_station("takaka", status = "open")
```

```{r, echo = FALSE}
takaka.df = structure(list(name = c("Takaka, Kotinga Road", "Riwaka At Takaka Hill", 
"Takaka Pohara", "Takaka At Harwoods", "Takaka At Kotinga", "Takaka @ Canaan", 
"Upper Takaka 2", "Takaka Ews", "Takaka Aero Raws"), network = c("F02882", 
"O12090", "F02884", "F15292", "F15291", "F0299A", "F12083", "F02885", 
"O00957"), agent = c(3788L, 44046L, 3790L, 44050L, 44051L, 44072L, 
11519L, 23849L, 41196L), start = structure(c(18273600, 316263600, 
520516800, 570020400, 704030400, 760014000, 805464000, 1020081600, 
1439294400), class = c("POSIXct", "POSIXt"), tzone = "NZ"), end = structure(c(1597665600, 
1597665600, 1597665600, 1597665600, 1597665600, 1597665600, 1597665600, 
1597665600, 1597665600), class = c("POSIXct", "POSIXt"), tzone = "NZ"), 
    open = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, 
    TRUE), distance = c(NA, NA, NA, NA, NA, NA, NA, NA, NA), 
    lat = c(-40.872, -41.03192, -40.845, -41.03094, -40.87068, 
    -40.93987, -41.01516, -40.86364, -40.81531), lon = c(172.809, 
    172.84439, 172.867, 172.79802, 172.808, 172.90821, 172.82582, 
    172.80568, 172.7765)), class = "data.frame", row.names = c(NA, 
-9L))

new("cfStation", takaka.df)
```

Stations are considered open in `clifro` if the final date returned from the
search is within four weeks of the current date. This gives the user a better
idea on the stations that are currently collecting data. 

## Station network ID search
The same can be done for searching stations using network ID although 
`search = "network"` needs to be added to the function call. Assume we knew
that the only stations we were interested in were the open stations whose 
network ID's match `F028`.

```{r, eval = FALSE}
cf_find_station("f028", search = "network", status = "all")
```

```{r, echo = FALSE}
xx.df = structure(list(name = c("Takaka, Kotinga Road", "Takaka Pohara", 
"Takaka Ews", "Aorere At Salisbury Bridge", "Takaka, Kotinga 2", 
"Nelson,Mckay Hut", "Gouland Downs", "Golden Bay,Table Hl I", 
"Golden Bay,Table Hl 2", "Tarakohe", "Takaka Aero", "Totaranui", 
"Takaka,Bu Bu", "Takaka", "Quartz Ranges"), network = c("F02882", 
"F02884", "F02885", "F02854", "F02883", "F02821", "F02831", "F02852", 
"F02853", "F02891", "F02871", "F02892", "F02872", "F02881", "F02851"
), agent = c(3788L, 3790L, 23849L, 44020L, 3789L, 3780L, 3781L, 
3783L, 3784L, 3791L, 3785L, 3792L, 3786L, 3787L, 3782L), start = structure(c(18273600, 
520516800, 1020081600, 1311595200, 502110000, 417960000, 467982000, 
233928000, 233928000, -1188819000, -1046692800, -410270400, -1159875000, 
-2082886200, -2177494200), class = c("POSIXct", "POSIXt"), tzone = "NZ"), 
    end = structure(c(1597665600, 1597665600, 1597665600, 1597665600, 
    1341057600, 745416000, 745416000, 690807600, 690807600, 599569200, 
    7732800, -294667200, -760190400, -1333452600, -2125049400
    ), class = c("POSIXct", "POSIXt"), tzone = "NZ"), open = c(TRUE, 
    TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
    FALSE, FALSE, FALSE, FALSE, FALSE), distance = c(NA, NA, 
    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA), lat = c(-40.872, 
    -40.845, -40.86364, -40.80236, -40.882, -40.89, -40.892, 
    -40.807, -40.807, -40.825, -40.816, -40.823, -40.85, -40.817, 
    -40.867), lon = c(172.809, 172.867, 172.80568, 172.53328, 
    172.801, 172.213, 172.351, 172.556, 172.556, 172.898, 172.772, 
    173.002, 172.733, 172.8, 172.517)), class = "data.frame", row.names = c(NA, 
-15L))

new("cfStation", xx.df)
```

Notice that the resulting dataframes in all of these searches are first ordered
by the date they last received data, and then by the date they opened, to return the
longest-running open stations first and the most historic, closed stations last.

## Return all stations within a region

This broad search returns all, open or closed stations within one of the 29 
preselected New Zealand regions (note that stations can belong to more than
one region). The `search = "region"` argument must be 
added to the `cf_find_station` function to conduct these searches. If the region 
is unknown then the search argument may be missing which brings up an 
interactive menu of the 29 regions for the user to select 
(`cf_find_station(search = "region")`), otherwise partial matching is used.

```{r, echo = FALSE}
open.queenstown.stations.df = dget(system.file("extdata", "queenStations", package = "clifro"))
open.queenstown.stations = new("cfStation", open.queenstown.stations.df)
```

```{r, eval = FALSE}
# Partial match for the Queenstown region
open.queenstown.stations = cf_find_station("queen", search = "region")
```

Typing `open.queenstown.stations` into R will then return all the 
`r nrow(open.queenstown.stations)` open Queenstown stations. This 
is clearly a burden to choose stations based on a large list of numbers hence 
plotting them on a map (covered below) to assess their spatial extent will make 
this task much easier.

## Return all stations within the vicinity of a given location

This location based search is conducted by including the 
`search = "latlong"` argument to the `cf_find_station` function. There are 
three parameters needed for this search; latitude, longitude and radius 
(kilometres). Just like any other function in R, if these arguments aren't 
named then the order matters and should be written in the order specified above.
The latitude and longitude must be given in decimal degrees.

We are (still) interested in finding all open stations around the small town of
Takaka. From 
[GeoHack](https://tools.wmflabs.org/geohack/geohack.php?pagename=Takaka%2C_New_Zealand&params=40_51_S_172_48_E_type:city%281149%29_region:NZ)
we can see that the latitude is -40.85 and the longitude is 172.8. We are 
interested in all open stations within a 10km radius of the main township.

```{r, echo = FALSE}
takaka.town.df = structure(list(name = c("Takaka, Kotinga Road", "Takaka Pohara", 
"Anatoki At Happy Sams", "Takaka At Kotinga", "Takaka Ews", "Motupiko At Reillys Bridge", 
"Takaka Aero Raws"), network = c("F02882", "F02884", "F15293", 
"F15291", "F02885", "F1529M", "O00957"), agent = c(3788L, 3790L, 
44015L, 44051L, 23849L, 44041L, 41196L), start = structure(c(18273600, 
520516800, 657284400, 704030400, 1020081600, 1164711600, 1439294400
), class = c("POSIXct", "POSIXt"), tzone = "NZ"), end = structure(c(1598788800, 
1598788800, 1598788800, 1598788800, 1598788800, 1598788800, 1598788800
), class = c("POSIXct", "POSIXt"), tzone = "NZ"), open = c(TRUE, 
TRUE, TRUE, TRUE, TRUE, TRUE, TRUE), distance = c(2.6, 5.7, 5.8, 
2.4, 1.6, 2.7, 4.3), lat = c(-40.872, -40.845, -40.88587, -40.87068, 
-40.86364, -40.85607, -40.81531), lon = c(172.809, 172.867, 172.74982, 
172.808, 172.80568, 172.83162, 172.7765)), class = "data.frame", row.names = c(NA, 
-7L))
takaka.town.st = new("cfStation", takaka.town.df)
```

```{r, eval = FALSE}
takaka.town.st = cf_find_station(lat = -40.85, long = 172.8, rad = 10, search = "latlong")

# Print the result, but remove the lat and lon columns to fit the page
takaka.town.st[, -c(8, 9)]
```

```{r, echo = -1}
takaka.town.st[, -c(8, 9)]

# We may rather order the stations by distance from the township
takaka.town.st[order(takaka.town.st$distance), -c(8, 9)]
```

# Searches based on datatypes

All the above searches did not include a datatype therefore they ignore the 
datatypes available at these stations. Imagine we are looking for 
hourly rain data at an open station in Takaka (using any of the aforementioned
searches), we would need to include the hourly rain datatype in the search for 
it to return a suitable station.

### Note
Unless the Reefton EWS station is the only CliFlo station of interest, the user 
will need a [CliFlo account](https://cliflo.niwa.co.nz/pls/niwp/wsubform.intro)
to get data from other stations.

```{r, echo = FALSE}
hourly.rain.dt = new("cfDatatype"
    , dt_name = "Precipitation"
    , dt_type = "Rain (fixed periods)"
    , dt_sel_option_names = list("Hourly")
    , dt_sel_combo_name = NA_character_
    , dt_param = structure("ls_ra,1,2,3,4", .Names = "dt1")
    , dt_sel_option_params = list(structure("182", .Names = "prm2"))
    , dt_selected_options = list(2)
    , dt_option_length = 4
)
```

```{r, eval = FALSE}
# Create a clifro datatype for hourly rain
hourly.rain.dt = cf_datatype(3, 1, 2)
hourly.rain.dt
```

```{r, echo = FALSE}
hourly.rain.dt
```

```{r, eval = FALSE}
# Conduct the search
cf_find_station("takaka", datatype = hourly.rain.dt)
```

```
##          name network agent      start        end open distance
## 1) Takaka Ews  F02885 23849 2002-06-02 2020-08-16 TRUE       NA
```

This tells us that the only *open* station in Takaka where hourly rain data 
is available is at the Takaka Ews station. 

# More than one search at a time

Since the `cf_find_station` function returns `cfStation` objects, any of these 
methods work on objects created from the `cf_station` function (see the 
[working with clifro stations vignette][clifrostation] for more details). We can 
conduct two or more searches at a time using the addition sign, just like we did 
for `cfDatatype`s (see the [choose datatypes vignette][chooseDatatype]).

We would like to return all open stations within a 10km radius of the Takaka 
township in the South Island, and the open stations in Kaitaia, in the North 
Island that collect hourly rain data.

```{r, echo = FALSE}
kaitaia.df = structure(list(name = c("Kaitaia Aero Ews", "Trounson Cws", "Russell Cws", 
"Kaikohe Aws", "Purerua Aws", "Cape Reinga Aws", "Kerikeri Aerodrome Aws", 
"Kaitaia Ews", "Dargaville 2 Ews", "Kerikeri Ews"), network = c("A53026", 
"A53762", "A54212", "A53487", "A54101", "A42462", "A53295", "A53127", 
"A53987", "A53191"), agent = c(18183L, 37131L, 41262L, 1134L, 
1196L, 1002L, 37258L, 17067L, 25119L, 1056L), start = structure(c(960984000, 
1244030400, 1459771200, 500727600, 788871600, 788871600, 1214395200, 
913806000, 1067425200, 1025179200), class = c("POSIXct", "POSIXt"
), tzone = "NZ"), end = structure(c(1598702400, 1598702400, 1598702400, 
1598616000, 1598616000, 1598616000, 1598616000, 1598443200, 1598011200, 
1597924800), class = c("POSIXct", "POSIXt"), tzone = "NZ"), open = c(TRUE, 
TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE), distance = c(NA, 
NA, NA, NA, NA, NA, NA, NA, NA, NA), lat = c(-35.0677, -35.72035, 
-35.26835, -35.4172, -35.129, -34.42963, -35.262, -35.13352, 
-35.93145, -35.183), lon = c(173.2874, 173.65153, 174.136, 173.8229, 
174.015, 172.68186, 173.911, 173.26294, 173.85317, 173.926)), class = "data.frame", row.names = c(NA, 
-10L))
kaitaia.st = new("cfStation", kaitaia.df)
my.composite.search = takaka.town.st + kaitaia.st
```

```{r, eval = FALSE}
my.composite.search = takaka.town.st + cf_find_station("kaitaia", 
                                                       search = "region", 
                                                       datatype = hourly.rain.dt)
my.composite.search
```

```{r, echo = -1}
my.composite.search

# How long have these stations been open for?
transform(my.composite.search, ndays = round(end - start))[, c(1, 10)]
```

# So where are these stations?

Up until now there probably hasn't been any good reason to choose clifro to 
search for stations instead of the 
['Choose Stations' form on CliFlo](https://cliflo.niwa.co.nz/pls/niwp/wstn.get_stn_html). 
However, the real advantage of using clifro is to visualise the station 
locations on a map by returning a KML file, particularly when there are lots of 
stations returned by the search. This Keyhole Markup Language 
([KML](https://resources.arcgis.com/en/help/main/10.1/index.html#//00s20000000m000000)) 
is an XML-based language provided by Google(TM) for defining the graphic display 
of spatial data in applications such as Google Earth(TM) and Google Maps(TM).

To return the stations as a KML file simply use the `cf_save_kml` function on 
any `cfStation` object. The `cf_find_station` function returns `cfStation` 
objects therefore it's very easy to plot these on a map. To assess the 
geographic extent of the Auckland stations we can return a KML file from the 
search and open it using our preferred KML-friendly software.

```{r,eval = FALSE}
# First, search for the stations
all.auckland.st = cf_find_station("auckland", search = "region", status = "all")
```

Now `all.auckland.st` contains the hundreds of Auckland stations where data have been recorded on CliFlo. 

```{r,eval=FALSE}
# Then save these as a KML
cf_save_kml(all.auckland.st, file_name = "all_auckland_stations")
```

The green markers represent the open stations and the red markers indicate 
closed stations. The resulting KML file is saved to the current R session's 
working directory by default. Have a look at the 
[clifro station vignette][clifrostation] for more methods and plotting of 
`cfStation` objects.

![All Auckland Stations][allAucklandStations]

[chooseDatatype]: choose-datatype.html
[clifrostation]: cfStation.html
[allAucklandStations]: figures/map.png
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findStations.R
\name{cf_save_kml}
\alias{cf_save_kml}
\title{Save Clifro Station Information to a KML File}
\usage{
cf_save_kml(station, file_name = "my_stations_", file_path = ".")
}
\arguments{
\item{station}{\code{cfStation} object containing one or more stations}

\item{file_name}{file name for the resulting KML file}

\item{file_path}{file path for the resulting KML file}
}
\description{
Save \code{\link{cfStation}} object information to a KML file.
}
\details{
The \code{cf_save_kml} function is for \code{\link{cfStation}}
objects to allow for the spatial visualisation of the selected stations. The
resulting KML file is saved and can then be opened by programs like Google
Earth (TM). The resultant KML file has the station names and locations shown
with green markers for open and red markers for closed stations. The agent
numbers, network ID's and date ranges are contained within the descriptions
for each station.

If no file name is specified, unique names are produced in the current \R
working directory.
}
\note{
The \code{.kml} suffix is appended automatically if it isn't already
present in the \code{file_name} argument.
}
\examples{
\dontrun{
# A selection of four Auckland region stations down the East Coast to the
# upper Waitemata Harbour; Leigh 2 Ews, Warkworth Ews, Tiri Tiri Lighthouse
# and Henderson
my.stations = cf_station(17838, 1340, 1401, 12327)
my.stations

# Save these stations to a KML file
cf_save_kml(my.stations)

# Double click on the file to open with a default program (if available). All
# the markers are green, indicating all these stations are open.

# Where is the subscription-free Reefton Ews station?
cf_save_kml(cf_station(), file_name = "reeftonEWS")

# It's located in the sou'west quadrant of Reefton town, in the upper, western
# part of the South Island, NZ.

# Find all the open and closed Christchurch stations (using partial matching)
all.chch.st = cf_find_station("christ", status = "all", search = "region")

# How many stations in total?
nrow(all.chch.st)

# Save all the Christchurch stations
cf_save_kml(all.chch.st, file_name = "all_Chch_stations")
}
}
\seealso{
\code{\link{cf_station}} and \code{vignette("cfStation")} for
working with stations when the agent numbers are known, otherwise
\code{\link{cf_find_station}} and code{vignette("choose-station")} for
creating \code{cfStation} objects when the agent numbers are unknown.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cfQuery.R
\name{cf_last_query}
\alias{cf_last_query}
\title{Retrieve Last Query Result from CliFlo}
\usage{
cf_last_query()
}
\description{
Retrieve the last query submitted to CliFlo instead of querying the database
again and losing subscription rows.
}
\details{
This function is a back up for when the clifro query has been submitted and
the data returned but has not been assigned, or inadvertently deleted. This
saves the user resubmitting queries and using more rows from their
subscription than needed.
}
\note{
Only the data from the last query is saved in \code{clifro}.
}
\examples{
\dontrun{
# Query CliFlo for wind at Reefton Ews
cf_query(cf_user(), cf_datatype(2, 1, 1, 1), cf_station(), "2012-01-01 00")

# Oops! Forgot to assign it to a variable...
reefton.wind = cf_last_query()
reefton.wind
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cfData-plotMethods.R
\name{windrose}
\alias{windrose}
\title{Plot a windrose}
\usage{
windrose(
  speed,
  direction,
  facet,
  n_directions = 12,
  n_speeds = 5,
  speed_cuts = NA,
  col_pal = "GnBu",
  ggtheme = c("grey", "gray", "bw", "linedraw", "light", "minimal", "classic"),
  legend_title = "Wind Speed",
  calm_wind = 0,
  variable_wind = 990,
  n_col = 1,
  ...
)
}
\arguments{
\item{speed}{numeric vector of wind speeds.}

\item{direction}{numeric vector of wind directions.}

\item{facet}{character or factor vector of the facets used to plot the various
windroses.}

\item{n_directions}{the number of direction bins to plot (petals on the rose).
The number of directions defaults to 12.}

\item{n_speeds}{the number of equally spaced wind speed bins to plot. This is
used if \code{speed_cuts} is \code{NA} (default 5).}

\item{speed_cuts}{numeric vector containing the cut points for the wind speed
intervals, or \code{NA} (default).}

\item{col_pal}{character string indicating the name of the
\code{RColorBrewer} colour palette to be
used for plotting, see 'Theme Selection' below.}

\item{ggtheme}{character string (partially) matching the
\code{\link[ggplot2]{ggtheme}} to be used for plotting, see
'Theme Selection' below.}

\item{legend_title}{character string to be used for the legend title.}

\item{calm_wind}{the upper limit for wind speed that is considered calm
(default 0).}

\item{variable_wind}{numeric code for variable winds (if applicable).}

\item{n_col}{The number of columns of plots (default 1).}

\item{...}{further arguments passed to \code{\link[ggplot2]{theme}}.}
}
\value{
a \code{ggplot} object.
}
\description{
Plot a windrose showing the wind speed and direction for given facets using
\pkg{ggplot2}.
}
\details{
This is intended to be used as a stand-alone function for any wind dataset. A
different windrose is plotted for each level of the faceting variable which
is coerced to a factor if necessary. The facets will generally be the station
where the data were collected, seasons or dates. Currently only one faceting
variable is allowed and is passed to \code{\link[ggplot2]{facet_wrap}} with
the formula \code{~facet}.
}
\section{Theme Selection}{

For black and white windroses that may be preferred if plots are to be used
in journal articles for example, recommended \code{ggtheme}s are \code{'bw'},
\code{'linedraw'}, \code{'minimal'} or \code{'classic'} and
the \code{col_pal} should be \code{'Greys'}. Otherwise, any of the sequential
\code{RColorBrewer} colour palettes are recommended for
colour plots.
}

\examples{
# Create some dummy wind data with predominant south to westerly winds, and
# occasional yet higher wind speeds from the NE (not too dissimilar to
# Auckland).

wind_df = data.frame(wind_speeds = c(rweibull(80, 2, 4), rweibull(20, 3, 9)),
                     wind_dirs = c(rnorm(80, 135, 55), rnorm(20, 315, 35)) \%\% 360,
                     station = rep(rep(c("Station A", "Station B"), 2),
                                   rep(c(40, 10), each = 2)))

# Plot a simple windrose using all the defaults, ignoring any facet variable
with(wind_df, windrose(wind_speeds, wind_dirs))

# Create custom speed bins, add a legend title, and change to a B&W theme
with(wind_df, windrose(wind_speeds, wind_dirs,
                       speed_cuts = c(3, 6, 9, 12),
                       legend_title = "Wind Speed\n(m/s)",
                       legend.title.align = .5,
                       ggtheme = "bw",
                       col_pal = "Greys"))

# Note that underscore-separated arguments come from the windrose method, and
# period-separated arguments come from ggplot2::theme().

# Include a facet variable with one level
with(wind_df, windrose(wind_speeds, wind_dirs, "Artificial Auckland Wind"))

# Plot a windrose for each level of the facet variable (each station)
with(wind_df, windrose(wind_speeds, wind_dirs, station, n_col = 2))

\dontrun{
# Save the plot as a png to the current working directory
library(ggplot2)
ggsave("my_windrose.png")
}

}
\seealso{
\code{\link[ggplot2]{theme}} for more possible arguments to pass to
\code{windrose}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cfData-plotMethods.R
\name{plot,cfRain,missing-method}
\alias{plot,cfRain,missing-method}
\alias{plot.cfRain}
\title{Plot Rain Time series}
\usage{
\S4method{plot}{cfRain,missing}(
  x,
  y,
  include_runoff = TRUE,
  ggtheme = c("grey", "gray", "bw", "linedraw", "light", "minimal", "classic"),
  scales = c("fixed", "free_x", "free_y", "free"),
  n_col = 1,
  ...
)
}
\arguments{
\item{x}{a \code{cfRain} object.}

\item{y}{missing.}

\item{include_runoff}{a logical indicating whether to plot the soil moisture
deficit and runoff as well as the rainfall, if the data
is available (default \code{TRUE}).}

\item{ggtheme}{character string (partially) matching the
\code{\link[ggplot2]{ggtheme}} to be used for plotting, see
'Theme Selection' below.}

\item{scales}{character string partially matching the \code{scales} argument
in the \code{link[ggplot2]{facet_wrap}} function.}

\item{n_col}{the number of columns of plots (default 1).}

\item{...}{further arguments passed to \code{\link[ggplot2]{theme}}.}
}
\description{
Plot the amount of rainfall (mm) through time, with optional available soil
water capacity and runoff amounts (if applicable).
}
\details{
When there is a rain event, the amount of runoff, if any, is dependent on how
much capacity the soil has available for more water. If there is no available
water capacity left in the soil then more rain will lead to a runoff event.
If \code{include_runoff = TRUE}, the available water capacity is plotted as
negative values and the runoff as positive values to signify this negative
relationship.
}
\examples{
\dontrun{
# Retrieve public rain data for a month from CliFlo (at Reefton Ews station)
reefton_rain = cf_query(cf_user(), cf_datatype(3, 1, 1), cf_station(),
                        start_date = "2012-08-01-00",
                        end_date = "2012-09-01-00")

class(reefton_rain) # cfRain object

# Plot the rain data using the defaults
plot(reefton_rain)

# Change the ggtheme and enlarge the text
library(ggplot2) # for element_text()
plot(reefton_rain, ggtheme = "bw", text = element_text(size = 16))

# Save the plot as a png to the current working directory
library(ggplot2) # for ggsave()
ggsave("my_rain_plot.png")
}
}
\seealso{
\code{\link{plot,cfDataList,missing-method}} for general
  information on default plotting of \code{cfData} and \code{cfDataList}
  objects, and the links within. See \code{\link{cf_query}} for creating
  \code{cfRain} objects.

  Refer to \code{\link[ggplot2]{theme}} for more possible arguments to pass
  to these methods.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cfDatatype.R
\name{cfDatatype-class}
\alias{cfDatatype-class}
\alias{cf_datatype}
\alias{cfDatatype}
\title{The Clifro Datatype Object}
\usage{
cf_datatype(
  select_1 = NA,
  select_2 = NA,
  check_box = NA,
  combo_box = NA,
  graphics = FALSE
)
}
\arguments{
\item{select_1}{a numeric vector of first node selections}

\item{select_2}{a numeric vector of second node selections}

\item{check_box}{a list containing the check box selections}

\item{combo_box}{a numeric vector containing the combo box selection
(if applicable)}

\item{graphics}{a logical indicating whether a graphics menu should be used,
if available}
}
\value{
\code{cfDatatype} object
}
\description{
Create a \code{cfDatatype} object by selecting one or more CliFlo datatypes
to build the \pkg{clifro} query.
}
\details{
An object inheriting from the \code{\link{cfDatatype}} class is created by
the constructor function \code{\link{cf_datatype}}. The function allows the
user to choose datatype(s) interactively (if no arguments are given), or to
create datatypes programmatically if the tree menu nodes are known a priori
(see examples). This function uses the same nodes, check box and combo box
options as CliFlo and can be viewed at the
\href{https://cliflo.niwa.co.nz/pls/niwp/wgenf.choose_datatype?cat=cat1}{datatype selection page}.
}
\note{
For the 'public' user (see examples) only the Reefton Ews station data
is available.

Currently clifro does not support datatypes from the special datasets
(Ten minute, Tier2, Virtual Climate, Lysimeter) or upper air measurements
from radiosondes and wind radar.
}
\examples{
\dontrun{
# Select the surface wind datatype manually (unknown tree nodes)
hourly.wind.dt = cf_datatype()
#  2  --> Datatype:        Wind
#  1  --> Datatype 2:      Surface Wind
#  2  --> Options:         Hourly Wind
# (2) --> Another option:  No
#  3  --> Units:           Knots
hourly.wind.dt

# Or select the datatype programatically (using the selections seen above)
hourly.wind.dt = cf_datatype(2, 1, 2, 3)
hourly.wind.dt
}
}
\seealso{
\code{\link{cf_user}} to create a \pkg{clifro} user,
  \code{\link{cf_station}} to choose the CliFlo stations and
  \code{vignette("choose-datatype")} for help choosing \code{cfDatatype}s.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cfUser.R
\docType{package}
\name{clifro}
\alias{clifro}
\alias{clifro-package}
\title{From CliFlo to \pkg{clifro}: Enhancing The National Climate Database With \R}
\description{
Import data from New Zealand's National Climate Database via CliFlo into \R
for exploring, analysis, plotting, exporting to KML, CSV, or other software.
}
\details{
The \pkg{clifro} package is intended to simplify the process of data
extraction, formatting and visualisation from the
\href{https://cliflo.niwa.co.nz/}{CliFlo web portal}. It
requires the user to build a query consisting of 3 main components; the user,
the datatype(s) and the station(s). These are
then combined using the \code{\link{cf_query}} function that sends the query
to the CliFlo database and returns the results that can easily be plotted
using generic plotting functions.

This package requires the user to already have a current subscription to the
National Climate Database unless a public user is sought, where data is
limited to Reefton Ews. Subscription is free and can obtained from
\url{https://cliflo.niwa.co.nz/pls/niwp/wsubform.intro}.
}
\examples{
\dontrun{
# Create a public user ----------------------------------------------------

public.user = cf_user() # Defaults to "public"
public.user

# Select datatypes --------------------------------------------------------

# 9am Surface wind (m/s)
wind.dt = cf_datatype(2, 1, 4, 1)

# Daily Rain
rain.dt = cf_datatype(3, 1, 1)

# Daily temperature extremes
temp.dt = cf_datatype(4, 2, 2)

# Combine them together
all.dts = wind.dt + rain.dt + temp.dt
all.dts

# Select the Reefton Ews station ------------------------------------------

reefton.st = cf_station()
reefton.st

# Submit the query --------------------------------------------------------

# Retrieve all data from ~ six months ago at 9am
reefton.data = cf_query(public.user, all.dts, reefton.st,
                        paste(as.Date(Sys.time()) - 182, "9"))
reefton.data


# Plot the data -----------------------------------------------------------

# Plot the 9am surface wind data (first dataframe in the list) ---
reefton.data[1]

# all identical - although passed to different methods
plot(reefton.data)    #plot,cfDataList,missing-method
plot(reefton.data, 1) #plot,cfDataList,numeric-method
plot(reefton.data[1]) #plot,cfData,missing-method --> plot,cfWind,missing-method

speed_plot(reefton.data)
direction_plot(reefton.data)

# Plot the daily rain data (second dataframe in the list) ---
reefton.data[2]

# With runoff and soil deficit
plot(reefton.data, 2)

# Just plot amount of rain (mm)
plot(reefton.data, 2, include_runoff = FALSE)

# Plot the hourly temperature data (third dataframe in the list) ---
plot(reefton.data, 3)

# Pass an argument to ggplot2::theme
library(ggplot2) # for element_text()
plot(reefton.data, 3, text = element_text(size = 18))
}
}
\seealso{
\code{\link{cf_user}}, \code{\link{cf_datatype}}, and
  \code{\link{cf_station}} for choosing the clifro user, datatypes and
  stations, respectively.
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cfUser.R
\name{cf_curl_opts}
\alias{cf_curl_opts}
\title{Store curl options for use within \pkg{clifro}}
\usage{
cf_curl_opts(..., .opts = list())
}
\arguments{
\item{...}{a name-value pairs that are passed to \code{RCurl curlOptions}}

\item{.opts}{a named list or \code{CURLOptions} object that are passed to \code{RCurl curlOptions}}
}
\description{
The \code{cf_curl_opts} function stores specific curl options that are used
for all the \pkg{clifro} queries.
}
\examples{
\dontrun{
# Specify options for use in all the curl handles created in clifro
cf_curl_opts(.opts = list(proxy = "http://xxxxx.yyyy.govt.nz:8080",
                          proxyusername  = "uid",
                          proxypassword  = "pwd",
                          ssl.verifypeer = FALSE))
# Or alternatively:
cf_curl_opts(proxy = "http://xxxxx.yyyy.govt.nz:8080",
             proxyusername  = "uid",
             proxypassword  = "pwd",
             ssl.verifypeer = FALSE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cfData-plotMethods.R
\name{summary,cfWind-method}
\alias{summary,cfWind-method}
\title{Summarise Clifro Wind Data}
\usage{
\S4method{summary}{cfWind}(object, calm_wind = 0)
}
\arguments{
\item{object}{a \code{cfWind} object.}

\item{calm_wind}{a single number containing the wind speed that is considered
calm.}
}
\description{
This is a summary method for \code{cfWind} objects.
}
\details{
A dataframe is returned containing the percentage of calm days
(wind speed >= \code{calm_days}), percentage of variable days (wind speed =
990), and quantiles from the empirical cumulative distribution functions for
each CliFlo station at which there is wind data.
}
\examples{
\dontrun{
# Retrieve maximum wind gust data at the Reefton Ews station from CliFlo
# (public data)
reefton_wind = cf_query(cf_user(), cf_datatype(2, 2, 1, 1), cf_station(),
                        start_date = "2012-01-01-00")

class(reefton_wind) # cfWind object

# Summarise the information
summary(reefton_wind)
}
}
\seealso{
\code{\link{plot.cfWind}} for default plotting of
 clifro wind data, and \code{\link{cf_query}} for creating \code{cfWind}
 objects.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cfQuery.R
\name{cf_query}
\alias{cf_query}
\title{Retrieve Data from the National Climate Database}
\usage{
cf_query(
  user,
  datatype,
  station,
  start_date,
  end_date = now(tz),
  date_format = "ymd_h",
  tz = "Pacific/Auckland",
  output_tz = c("local", "NZST", "UTC"),
  quiet = FALSE
)
}
\arguments{
\item{user}{a \code{\link{cfUser}} object.}

\item{datatype}{a \code{\link{cfDatatype}} object containing the datatypes to
be retrieved.}

\item{station}{a \code{\link{cfStation}} object containing the stations where
the datatypes will be retrieved from.}

\item{start_date}{a character, Date or POSIXt object indicating the start
date. If a character string is supplied the date format
should be in the form \code{yyyy-mm-dd-hh} unless
\code{date_format} is specified.}

\item{end_date}{same as \code{start_date}. Defaults to
\code{\link[lubridate]{now}}.}

\item{date_format}{a character string matching one of \code{"ymd_h"},
\code{"mdy_h"}, \code{"ydm_h"} or \code{"dmy_h"}
representing the \code{\link[lubridate]{lubridate-package}}
date parsing function.}

\item{tz}{the timezone for which the start and end dates refer to. Conversion
to Pacific/Auckland time is done automatically through the
\code{\link[lubridate]{with_tz}} function. Defaults to
"Pacific/Auckland".}

\item{output_tz}{the timezone of the output. This can be one of either "local",
"UTC", or "NZST".}

\item{quiet}{logical. When \code{TRUE} the function evaluates without
displaying customary messages. Messages from CliFlo are still
displayed.}
}
\value{
a \code{cfData} or \code{cfDataList} object.
}
\description{
Query the National Climate Database via CliFlo based on the \pkg{clifro} user
and selected datatypes, stations and dates.
}
\details{
The \code{cf_query} function is used to combine the \pkg{clifro} user
(\code{\link{cfUser}}), along with the desired datatypes
(\code{\link{cfDatatype}}) and stations (\code{\link{cfStation}}). The query
is 'built up' using these objects, along with the necessary dates. The
function then uses all these data to query the National Climate Database via
the CliFlo web portal and returns one of the many \code{cfData}
objects if one dataframe is returned, or a \code{cfDataList} object if
there is more than one dataframe returned from CliFlo. If a \code{cfDataList}
is returned, each element in the list is a subclass of the \code{cfData}
class, see the 'cfData Subclasses' section.
}
\section{CfData Subclasses}{


There are 8 \code{cfData} subclasses that are returned from \code{cf_query}
depending on the datatype requested. Each of these subclasses have default
\code{plot} methods for usability and efficiency in exploring and plotting
\pkg{clifro} data.

The following table summarises these subclasses and how they are created, see
also the examples on how to automatically create some of these subclasses.

\tabular{ll}{
\strong{Subclass} \tab \strong{CliFlo Datatype}\cr
cfWind \tab Any 'Wind' data \cr
cfRain \tab Any 'Precipitation' data \cr
cfScreen Obs \tab 'Temperature and Humidity' data measured in a standard screen \cr
cfTemp \tab Maximum and minimum 'Temperature and Humidity' data \cr
cfEarthTemp \tab 'Temperature and Humidity' data at a given depth \cr
cfSunshine \tab Any 'Sunshine & Radiation' data \cr
cfPressure \tab Any 'Pressure' data \cr
cfOther \tab Any other CliFlo 'Daily and Hourly Observations' \cr
}
}

\examples{
\dontrun{
# Retrieve daily rain data from Reefton Ews
daily.rain = cf_query(cf_user("public"), cf_datatype(3, 1, 1),
                      cf_station(), "2012-01-01 00")
daily.rain

# returns a cfData object as there is only one datatype
class(daily.rain) # 'cfRain' object - inherits 'cfData'

# Look up the help page for cfRain plot methods
?plot.cfRain

# Retrieve daily rain and wind data from Reefton Ews

daily.dts = cf_query(cf_user("public"),
                     cf_datatype(c(2, 3), c(1, 1), list(4, 1), c(1, NA)),
                     cf_station(), "2012-01-01 00", "2013-01-01 00")
daily.dts

# returns a cfDataList object as there is more than one datatype. Each
# element of the cfDataList is an object inheriting from the cfData class.
class(daily.dts)     # cfDataList
class(daily.dts[1])  # cfRain
class(daily.dts[2])  # cfWind

# Create a cfSunshine object (inherits cfData)
# Retrieve daily global radiation data at Reefton Ews
rad.data = cf_query(cf_user(), cf_datatype(5,2,1), cf_station(),
                    "2012-01-01 00")
rad.data

# The cf_query function automatically creates the appropriate cfData subclass
class(rad.data)

# The advantage of having these subclasses is that it makes plotting very easy
plot(rad.data)
plot(daily.rain)
plot(daily.rain, include_runoff = FALSE)
plot(daily.dts)
plot(daily.dts, 2)
}
}
\seealso{
\code{\link{cf_user}}, \code{\link{cf_datatype}} and
  \code{\link{cf_station}} for creating the objects needed for a query. See
  \code{\link{plot,cfDataList,missing-method}} for general information on
  default plotting of \code{cfData} and \code{cfDataList} objects, and the
  links within.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cfData-plotMethods.R
\name{plot,cfScreenObs,missing-method}
\alias{plot,cfScreenObs,missing-method}
\alias{plot.cfScreenObs}
\title{Plot Screen Observations}
\usage{
\S4method{plot}{cfScreenObs,missing}(
  x,
  y,
  ggtheme = c("grey", "gray", "bw", "linedraw", "light", "minimal", "classic"),
  scales = c("fixed", "free_x", "free_y", "free"),
  n_col = 1,
  ...
)
}
\arguments{
\item{x}{a cfScreenObs object.}

\item{y}{missing.}

\item{ggtheme}{character string (partially) matching the
\code{\link[ggplot2]{ggtheme}} to be used for plotting, see
'Theme Selection' below.}

\item{scales}{character string partially matching the \code{scales} argument
in the \code{link[ggplot2]{facet_wrap}} function.}

\item{n_col}{the number of columns of plots (default 1).}

\item{...}{further arguments passed to \code{\link[ggplot2]{theme}}.}
}
\description{
Plot temperature data from screen observations (degrees celsius) through time.
}
\details{
Temperature data from screen observations include the air, and wet bulb,
temperature at the time the measurement was taken (dry bulb and wet bulb
respectively), and the dew point. The dew point is the air temperature at
which dew starts to form. That is the temperature to which a given air parcel
must be cooled at constant pressure and constant water vapour content in
order for saturation to occur.

The resulting figure plots the dry bulb, wet bulb and dew point temperatures
on the same  scale, for each station.
}
\examples{
\dontrun{
# Retrieve public temperature data from screen observations for the last week
# at Reefton Ews station

# Subtract 7 days from today's date to get the start date
last_week = paste(as.character(Sys.Date() - 7), 0)

reefton_screenobs = cf_query(cf_user(), cf_datatype(4, 1, 1), cf_station(),
                             start_date = last_week)

class(reefton_screenobs) # cfScreenObs object

# Plot the temperature data using the defaults
plot(reefton_screenobs)

# Enlarge the text and add the observations as points
library(ggplot2) # for element_text() and geom_point()
plot(reefton_screenobs, ggtheme = "bw", text = element_text(size = 16)) +
  geom_point(size = 3, shape = 1)

# Save the plot as a png to the current working directory
library(ggplot2) # for ggsave()
ggsave("my_screenobs_plot.png")
}
}
\references{
\href{https://cliflo.niwa.co.nz/pls/niwp/wh.do_help?id=ls_scr1}{Screen Observation details}.
}
\seealso{
\code{\link{plot,cfDataList,missing-method}} for general
  information on default plotting of \code{cfData} and \code{cfDataList}
  objects, and the links within. See \code{\link{cf_query}} for creating
  \code{cfScreenObs} objects.

  Refer to \code{\link[ggplot2]{theme}} for more possible arguments to pass
  to these methods.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cfUser.R
\name{valid_cfuser}
\alias{valid_cfuser}
\alias{cf_login}
\alias{cf_logout}
\title{Validation Functions For The \code{cfUser} Class}
\usage{
cf_login(object, ...)

cf_logout(object, msg = TRUE, ...)

valid_cfuser(object)
}
\arguments{
\item{object}{S4 object which inherits the \code{cfUser} class}

\item{...}{Other options passed to the \code{\link[httr]{GET}} or \code{\link[httr]{POST}} functions.}

\item{msg}{Display a 'successful logout' message, defaults to
\code{TRUE}.}
}
\description{
These internal functions are used by the \code{\link{cf_user}} constructor
function to ensure the user has a valid subscription to CliFlo.
}
\details{
\code{cf_login} initiates a curl handle storing the cookies in the current
\R session's temporary directory. It then POSTs the user credentials to the
CliFlo login page and stores the resultant \code{h1} heading to check for the
string 'Info'. The cookies are kept for future (immediate) use.

\code{cf_logout} points the curl handle to the existing cookie session
initiated with \code{cf_login}. It reads the header information from the
cliflo logout page to ensure no HTTP error and logs the user out on
cliflo and deletes the cookies. This should be (is) called immediately after
\code{cf_login} in any function requiring a login, using
\code{\link{on.exit}} to ensure the user isn't still logged in on the server,
after the function call, for any reason.

\code{valid_cfuser} is the validation function for the \code{cfUser} class
and uses  \code{cf_login} to ensure the credentials are authenticated on the
CliFlo server and then (\code{cf_})logs out immediately afterwards. It also
ensures the user provides exactly one username and password - except for
'public' users.
}
\examples{
\dontrun{
cf_user("public")                    # Returns a valid object
cf_user("bad_name", "bad_password")    # Bad Login
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cfStation.R, R/cfDatatype.R
\name{+,cfStation,cfStation-method}
\alias{+,cfStation,cfStation-method}
\alias{+,cfDatatype,cfDatatype-method}
\title{Arithmetic Operators for Clifro Objects}
\usage{
\S4method{+}{cfStation,cfStation}(e1, e2)

\S4method{+}{cfDatatype,cfDatatype}(e1, e2)
}
\arguments{
\item{e1}{a \code{cfDatatype} or \code{cfStation} object}

\item{e2}{an object matching the class of e1}
}
\description{
This operator allows you to add more datatypes or stations to
\code{cfDatatype} and \code{cfStation} objects respectively.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cfData-plotMethods.R
\name{plot,cfPressure,missing-method}
\alias{plot,cfPressure,missing-method}
\alias{plot.cfPressure}
\title{Plot Mean Sea Level Atmospheric Pressure}
\usage{
\S4method{plot}{cfPressure,missing}(
  x,
  y,
  ggtheme = c("grey", "gray", "bw", "linedraw", "light", "minimal", "classic"),
  scales = c("fixed", "free_x", "free_y", "free"),
  n_col = 1,
  ...
)
}
\arguments{
\item{x}{a cfPressure object.}

\item{y}{missing.}

\item{ggtheme}{character string (partially) matching the
\code{\link[ggplot2]{ggtheme}} to be used for plotting, see
'Theme Selection' below.}

\item{scales}{character string partially matching the \code{scales} argument
in the \code{link[ggplot2]{facet_wrap}} function.}

\item{n_col}{the number of columns of plots (default 1).}

\item{...}{further arguments passed to \code{\link[ggplot2]{theme}}.}
}
\description{
Plot the MSL atmospheric pressure through time.
}
\examples{
\dontrun{
# Retrieve public hourly atmospheric pressure data for the last 30 days at
# Reefton Ews station

# Subtract 30 days from today's date to get the start date
last_month = paste(as.character(Sys.Date() - 30), 0)

reefton_pressure = cf_query(cf_user(), cf_datatype(7, 1, 1), cf_station(),
                            start_date = last_month)

class(reefton_pressure) # cfPressure object

# Plot the atmospheric pressure data using the defaults
plot(reefton_pressure)

# Enlarge the text and add the observations as points
library(ggplot2) # for element_text() and geom_point()
plot(reefton_pressure, ggtheme = "bw", text = element_text(size = 16)) +
  geom_point(size = 3, shape = 1)

# Save the plot as a png to the current working directory
library(ggplot2) # for ggsave()
ggsave("my_pressure_plot.png")
}
}
\seealso{
\code{\link{plot,cfDataList,missing-method}} for general
  information on default plotting of \code{cfData} and \code{cfDataList}
  objects, and the links within. See \code{\link{cf_query}} for creating
  \code{cfPressure} objects.

  Refer to \code{\link[ggplot2]{theme}} for more possible arguments to pass
  to these methods.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cfData-plotMethods.R
\name{plot,cfWind,missing-method}
\alias{plot,cfWind,missing-method}
\alias{plot.cfWind}
\alias{direction_plot,cfWind,missing-method}
\alias{direction_plot}
\alias{direction_plot,cfDataList,missing-method}
\alias{direction_plot,cfDataList,numeric-method}
\alias{speed_plot,cfWind,missing-method}
\alias{speed_plot}
\alias{speed_plot,cfDataList,missing-method}
\alias{speed_plot,cfDataList,numeric-method}
\title{Plot Clifro Wind Objects}
\usage{
\S4method{plot}{cfWind,missing}(
  x,
  y,
  n_directions = 12,
  n_speeds = 5,
  speed_cuts = NULL,
  col_pal = "GnBu",
  ggtheme = c("grey", "gray", "bw", "linedraw", "light", "minimal", "classic"),
  n_col = 1,
  ...
)

\S4method{direction_plot}{cfWind,missing}(
  x,
  y,
  ggtheme = c("grey", "gray", "bw", "linedraw", "light", "minimal", "classic"),
  contours = 10,
  n_col = 1,
  ...
)

\S4method{direction_plot}{cfDataList,missing}(x, y, ...)

\S4method{direction_plot}{cfDataList,numeric}(x, y, ...)

\S4method{speed_plot}{cfWind,missing}(
  x,
  y,
  ggtheme = c("grey", "gray", "bw", "linedraw", "light", "minimal", "classic"),
  scales = c("fixed", "free_x", "free_y", "free"),
  n_col = 1,
  ...
)

\S4method{speed_plot}{cfDataList,missing}(x, y, ...)

\S4method{speed_plot}{cfDataList,numeric}(x, y, ...)
}
\arguments{
\item{x}{a \code{cfWind} or \code{cfDataList} object.}

\item{y}{missing if \code{x} is a .\code{cfWind} object, otherwise a number
indicating the dataframe to plot in the \code{cfDataList} (defaults
to 1).}

\item{n_directions}{the number of direction bins to plot (petals on the
rose). The number of directions defaults to 12.}

\item{n_speeds}{the number of equally spaced wind speed bins to plot. This is
used if \code{spd_cuts} is NA (default 5).}

\item{speed_cuts}{numeric vector containing the cut points for the wind speed
intervals, or NA (default).}

\item{col_pal}{character string indicating the name of the
\code{RColorBrewer} colour palette to be
used for plotting, see 'Theme Selection' below.}

\item{ggtheme}{character string (partially) matching the
\code{\link[ggplot2]{ggtheme}} to be used for plotting, see
'Theme Selection' below.}

\item{n_col}{the number of columns of plots (default 1).}

\item{...}{further arguments passed to \code{\link[ggplot2]{theme}}.}

\item{contours}{the number of contour lines to draw (default 10).}

\item{scales}{character string partially matching the \code{scales} argument
in the \code{link[ggplot2]{facet_wrap}} function.}
}
\description{
Various plot methods for exploring wind speed and direction patterns for
given CliFlo stations.
}
\details{
If \code{x} is a \code{cfDataList}, by default the first datatype will be
plotted unless \code{y} is supplied.
}
\note{
If \code{x} is a \code{cfDataList} object and \code{y} refers to a
\pkg{clifro} dataframe that is not a \code{cfWind} object then it will be
passed to another method, if available.

The default \code{plot} method plots a different windrose for each CliFlo
station. The \code{direction_plot} method plots wind direction contours
through time to visualise temporal patterns in wind directions. The
\code{speed_plot} method plots the time series of wind speeds with a +/-
standard deviation region (if applicable).

Given a value on the x-axis, the ends of the density function along the
 y-axis are not constrained to be equal for any of the derivatives for the
 \code{direction_plot} method. That is, the contours at direction = 0, do not
 match the contours at direction = 360.

 @seealso \code{\link{plot,cfDataList,missing-method}} for general
  information on default plotting of \code{cfData} and \code{cfDataList}
  objects, and the links within. See \code{\link{cf_query}} for creating
  \code{cfWind} objects or \code{\link{windrose}} for plotting any wind data.
  Refer to \code{\link[ggplot2]{theme}} for more possible arguments to pass
  to these methods. \code{\link{summary,cfWind-method}} for summarising wind
  information at each CliFlo station.
}
\section{Theme Selection}{

For black and white windroses that may be preferred if plots are to be used
in journal articles for example, recommended \code{ggtheme}s are \code{'bw'},
\code{'linedraw'}, \code{'minimal'} or \code{'classic'} and
the \code{col_pal} should be \code{'Greys'}. Otherwise, any of the sequential
\code{RColorBrewer} colour palettes are recommended for
colour plots.
}

\examples{
\dontrun{
# Retrieve maximum wind gust data at the Reefton Ews station from CliFlo
# (public data)
reefton_wind = cf_query(cf_user(), cf_datatype(2, 2, 1, 1), cf_station(),
                        start_date = "2012-01-01-00")

class(reefton_wind)

# Examples of the default plots --------------------------------------------

# Plot a windrose
plot(reefton_wind)

# Plot the wind direction contours
direction_plot(reefton_wind)

# Plot the wind speed time-series
speed_plot(reefton_wind)

# Examples of changing the defaults ----------------------------------------

# Plot black and white windroses
plot(reefton_wind, ggtheme = "bw", col_pal = "Greys")
plot(reefton_wind, ggtheme = "linedraw", col_pal = "Greys")
plot(reefton_wind, ggtheme = "classic", col_pal = "Greys")
plot(reefton_wind, ggtheme = "minimal", col_pal = "Greys")

# Plot the wind directions using 20 contours and the ggtheme 'classic'
direction_plot(reefton_wind, ggtheme = "classic", contours = 20)

# Enlarge all the text to 18pt
library(ggplot2) # for element_text() and geom_point()
direction_plot(reefton_wind, ggtheme = "classic", contours = 20,
               text = element_text(size = 18))

# Include the actual observations in the plots
direction_plot(reefton_wind) + geom_point(alpha = .2, size = 3)

speed_plot(reefton_wind, ggtheme = "classic", text = element_text(size = 16)) +
  geom_point(shape = 1, size = 3)
# or equivalently using base graphics:
plot(reefton_wind$Date, reefton_wind$Speed, type = 'o',
     xlab = NA, ylab = "Daily max gust (m/s)", las = 1, main = "Reefton Ews")

# Example of plotting a cfDataList -----------------------------------------
# Collect both surface wind run and hourly surface wind observations from
# Reefton Ews
reefton_list = cf_query(cf_user(), cf_datatype(2, 1, 1:2, 1),
                        cf_station(), "2012-01-01 00", "2012-02-01 00")

reefton_list

class(reefton_list) #cfDataList

# Plot the first (default) dataframe
plot(reefton_list) # Error - no wind directions for wind run datatypes
# Try speed_plot instead
speed_plot(reefton_list)

# Plot the second dataframe in the cfDataList
plot(reefton_list, 2)           # identical to plot(reefton_list[2])
speed_plot(reefton_list, 2)     # identical to speed_plot(reefton_list[2])
direction_plot(reefton_list, 2) # identical to direction_plot(reefton_list[2])

# Save the ggplot externally -----------------------------------------------

# Save the plot as a png to the current working directory
library(ggplot2) # for ggsave()
ggsave("my_wind_plot.png")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataFrame.R, R/cfStation.R, R/cfDataList.R,
%   R/cfDatatype.R
\docType{methods}
\name{[[,dataFrame-method}
\alias{[[,dataFrame-method}
\alias{[,dataFrame,ANY,ANY,ANY-method}
\alias{$,dataFrame-method}
\alias{[,cfStation,ANY,ANY,ANY-method}
\alias{Extract}
\alias{[,cfDataList,ANY,ANY,ANY-method}
\alias{[[,cfDataList-method}
\alias{[,cfDatatype,ANY,missing,missing-method}
\alias{[,cfDatatype,ANY,missing,missing}
\title{Subsetting Methods for Clifro Objects}
\usage{
\S4method{[[}{dataFrame}(x, i)

\S4method{[}{dataFrame,ANY,ANY,ANY}(x, i, j, drop)

\S4method{$}{dataFrame}(x, name)

\S4method{[}{cfStation,ANY,ANY,ANY}(x, i, j, drop = TRUE)

\S4method{[}{cfDataList,ANY,ANY,ANY}(x, i, j)

\S4method{[[}{cfDataList}(x, i)

\S4method{[}{cfDatatype,ANY,missing,missing}(x, i, j, drop)
}
\arguments{
\item{x}{a \pkg{clifro} object}

\item{i}{indices specifying elements to extract. Indices are
\code{numeric} or \code{character} vectors or empty (missing) or
\code{NULL}. Character vectors will be matched to the names of
the object.}

\item{j}{indices specifying elements to extract. Indices are
\code{numeric} or \code{character} vectors or empty (missing) or
\code{NULL}. Character vectors will be matched to the names of
the object.}

\item{drop}{if \code{TRUE}, the result is coerced to the lowest possible
dimension. See \code{\link{drop}} for further details.}

\item{name}{a literal character string. This is partially matched to the
names of the object.}
}
\description{
Operators acting on \code{cfDataList}, \code{cfDatatype}, \code{cfStation},
and \code{dataFrame} objects.
}
\details{
These are methods for the generic operators for classes within \pkg{clifro}.
They are intended to give the user the familiar functionality of subsetting
\code{\link{data.frame}} objects.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cfData-plotMethods.R
\name{plot,cfEarthTemp,missing-method}
\alias{plot,cfEarthTemp,missing-method}
\alias{plot.cfEarthTemp}
\title{Plot Earth Temperatures}
\usage{
\S4method{plot}{cfEarthTemp,missing}(
  x,
  y,
  ggtheme = c("grey", "gray", "bw", "linedraw", "light", "minimal", "classic"),
  scales = c("fixed", "free_x", "free_y", "free"),
  n_col = 1,
  ...
)
}
\arguments{
\item{x}{a cfEarthTemp object.}

\item{y}{missing.}

\item{ggtheme}{character string (partially) matching the
\code{\link[ggplot2]{ggtheme}} to be used for plotting, see
'Theme Selection' below.}

\item{scales}{character string partially matching the \code{scales} argument
in the \code{link[ggplot2]{facet_wrap}} function.}

\item{n_col}{the number of columns of plots (default 1).}

\item{...}{further arguments passed to \code{\link[ggplot2]{theme}}.}
}
\description{
Plot the earth temperature for a given depth (degrees celsius) through time,
for each chosen CliFlo station.
}
\examples{
\dontrun{
# Retrieve public earth temperature data for the last 30 days at Reefton Ews
# station, at a depth of 10cm

# Subtract 30 days from today's date to get the start date
last_month = paste(as.character(Sys.Date() - 30), 0)

reefton_earth = cf_query(cf_user(), cf_datatype(4, 3, 2), cf_station(),
                         start_date = last_month)

class(reefton_earth) # cfTemp object

# Plot the temperature data using the defaults
plot(reefton_earth)

# Enlarge the text and add the observations as points
library(ggplot2) # for element_text() and geom_point()
plot(reefton_earth, ggtheme = "bw", text = element_text(size = 16)) +
  geom_point(size = 3, shape = 1)

# Save the plot as a png to the current working directory
library(ggplot2) # for ggsave()
ggsave("my_earthTemp_plot.png")
}
}
\seealso{
\code{\link{plot,cfDataList,missing-method}} for general
  information on default plotting of \code{cfData} and \code{cfDataList}
  objects, and the links within. See \code{\link{cf_query}} for creating
  \code{cfEarthTemp} objects.

  Refer to \code{\link[ggplot2]{theme}} for more possible arguments to pass
  to these methods.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataFrame.R
\name{dimnames,dataFrame-method}
\alias{dimnames,dataFrame-method}
\alias{dim,dataFrame-method}
\title{Dimension Attributes of a Clifro Object}
\usage{
\S4method{dimnames}{dataFrame}(x)

\S4method{dim}{dataFrame}(x)
}
\arguments{
\item{x}{a \code{dataFrame} object

Specifically, a \code{dataFrame} object is any \code{\link{cfStation}} or 
\code{cfData} object. These functions are provided for the user to have (some)
familiar \code{data.frame}-type functions available for use on \pkg{clifro}
objects.}
}
\description{
Retrieve the dimensions or dimension names of a \code{dataFrame} object.
}
\seealso{
\code{\link{cf_query}} for creating \code{cfData} objects, and
  \code{\link{cf_station}} for creating \code{cfStation} objects.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cfUser.R
\docType{class}
\name{cfUser-class}
\alias{cfUser-class}
\alias{cf_user}
\alias{cfUser}
\title{The Clifro User Object}
\usage{
cf_user(username = "public", password = character())
}
\arguments{
\item{username}{a character string to be used as the cliflo username}

\item{password}{a character string to be used as the cliflo password}
}
\value{
\code{cfUser} object
}
\description{
Create a \code{cfUser} object to allow the user to log into CliFlo from \R
and  build their query.
}
\details{
An object inheriting from the \code{cfUser} class is created by the constructor
function \code{cf_user}. The user must have an active subscription to cliflo
in order to create a valid object, unless a 'public' user is sought.
Visit \url{https://cliflo.niwa.co.nz/} for more information and to subscribe
to cliflo.
}
\note{
For the 'public' user (see examples) only the Reefton Ews station data
is available.
}
\examples{
\dontrun{
public.cfuser = cf_user(username = "public")
public.cfuser
}
}
\seealso{
\code{\link{valid_cfuser}} for details on the validation of
\code{cfUser} and \code{\link{summary,cfUser-method}} to summarise user
information.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cfStation.R
\name{cfStation-class}
\alias{cfStation-class}
\alias{cf_station}
\alias{cfStation}
\title{The Clifro Station Object}
\usage{
cf_station(...)
}
\arguments{
\item{...}{comma separated agent numbers}
}
\value{
\code{cfStation} object
}
\description{
Create a \code{cfStation} object containing station information for one or
more CliFlo stations.
}
\details{
A \code{cfStation} object is created by the constructor function
\code{cf_station}. The unique agent numbers of the stations are all that is
required to create a \code{cfStation} object using the \code{cf_station}
function. The rest of the station information including the name, network and
agent ID, start and end dates, coordinates, as well as other data is scraped
from CliFlo.

This function is used for when the agent numbers are already known. For help
creating \code{cfStation} objects when the agent numbers are unknown see the
\code{\link{cf_find_station}} function.
}
\examples{
\dontrun{
# Create a cfStation object for the Leigh 1 and 2 Ews stations
leigh.st = cf_station(1339, 1340)
leigh.st

# Note, this can also be achieved using the '+' operator
leigh.st = cf_station(1339) + cf_station(1340)
leigh.st

# Add another column showing how long the stations have been open for
leigh.df = as(leigh.st, "data.frame")
leigh.df$ndays = with(leigh.df, round(end - start))
leigh.df

# Save the stations to the current working directory as a KML to visualise
# the station locations
cf_save_kml(leigh.st)
}
}
\seealso{
\code{\link{cf_find_station}} for creating \code{cfStation} objects
when the agent numbers are not known and \code{vignette("cfStation")}
for working with clifro stations including spatial plotting in \R. For saving
\code{cfStation} objects as KML files refer to the vignette or
\code{\link{cf_save_kml}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cfData-plotMethods.R
\name{plot,cfTemp,missing-method}
\alias{plot,cfTemp,missing-method}
\alias{plot.cfTemp}
\title{Plot Temperature Range}
\usage{
\S4method{plot}{cfTemp,missing}(
  x,
  y,
  ggtheme = c("grey", "gray", "bw", "linedraw", "light", "minimal", "classic"),
  scales = c("fixed", "free_x", "free_y", "free"),
  n_col = 1,
  ...
)
}
\arguments{
\item{x}{a cfTemp object.}

\item{y}{missing.}

\item{ggtheme}{character string (partially) matching the
\code{\link[ggplot2]{ggtheme}} to be used for plotting, see
'Theme Selection' below.}

\item{scales}{character string partially matching the \code{scales} argument
in the \code{link[ggplot2]{facet_wrap}} function.}

\item{n_col}{the number of columns of plots (default 1).}

\item{...}{further arguments passed to \code{\link[ggplot2]{theme}}.}
}
\description{
Plot minimum and maximum temperature data for a given period (degrees
celsius) through time, for each chosen CliFlo station.
}
\details{
This plotting method shows the temperature extremes as a grey region on the
plot, with a black line indicating the average temperature (if available).
}
\examples{
\dontrun{
# Retrieve public hourly minimum and maximum temperature data for the last
week at Reefton Ews station

# Subtract 7 days from today's date to get the start date
last_week = paste(as.character(Sys.Date() - 7), 0)

reefton_temp = cf_query(cf_user(), cf_datatype(4, 2, 2), cf_station(),
                        start_date = last_week)

class(reefton_temp) # cfTemp object

# Plot the temperature data using the defaults
plot(reefton_temp)

# Enlarge the text and add the observations as points
library(ggplot2) # for element_text() and geom_point()
plot(reefton_temp, ggtheme = "bw", text = element_text(size = 16)) +
  geom_point(size = 3, shape = 1)

# Save the plot as a png to the current working directory
library(ggplot2) # for ggsave()
ggsave("my_temperature_plot.png")
}
}
\seealso{
\code{\link{plot,cfDataList,missing-method}} for general
  information on default plotting of \code{cfData} and \code{cfDataList}
  objects, and the links within. See \code{\link{cf_query}} for creating
  \code{cfTemp} objects.

  Refer to \code{\link[ggplot2]{theme}} for more possible arguments to pass
  to these methods.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cfData-plotMethods.R
\name{plot,cfSunshine,missing-method}
\alias{plot,cfSunshine,missing-method}
\alias{plot.cfSunshine}
\title{Plot Sunshine Hours}
\usage{
\S4method{plot}{cfSunshine,missing}(
  x,
  y,
  ggtheme = c("grey", "gray", "bw", "linedraw", "light", "minimal", "classic"),
  scales = c("fixed", "free_x", "free_y", "free"),
  n_col = 1,
  ...
)
}
\arguments{
\item{x}{a cfSunshine object.}

\item{y}{missing.}

\item{ggtheme}{character string (partially) matching the
\code{\link[ggplot2]{ggtheme}} to be used for plotting, see
'Theme Selection' below.}

\item{scales}{character string partially matching the \code{scales} argument
in the \code{link[ggplot2]{facet_wrap}} function.}

\item{n_col}{the number of columns of plots (default 1).}

\item{...}{further arguments passed to \code{\link[ggplot2]{theme}}.}
}
\description{
Plot the duration of accumulated bright sunshine hours through time.
}
\examples{
\dontrun{
# Retrieve public hourly sunshine data for the last 7 days at Reefton Ews
# station

# Subtract 7 days from today's date to get the start date
last_week = paste(as.character(Sys.Date() - 7), 0)

reefton_sun = cf_query(cf_user(), cf_datatype(5, 1, 2), cf_station(),
                       start_date = last_week)

class(reefton_sun) # cfSunshine object

# Plot the temperature data using the defaults
plot(reefton_sun)

# Enlarge the text and add the observations as points
library(ggplot2) # for element_text() and geom_point()
plot(reefton_sun, ggtheme = "bw", text = element_text(size = 16)) +
  geom_point(size = 3, shape = 1)

# Save the plot as a png to the current working directory
library(ggplot2) # for ggsave()
ggsave("my_sunshine_plot.png")
}
}
\seealso{
\code{\link{plot,cfDataList,missing-method}} for general
  information on default plotting of \code{cfData} and \code{cfDataList}
  objects, and the links within. See \code{\link{cf_query}} for creating
  \code{cfSunshine} objects.

  Refer to \code{\link[ggplot2]{theme}} for more possible arguments to pass
  to these methods.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cfData-plotMethods.R
\name{plot.cfDataList}
\alias{plot.cfDataList}
\alias{plot,cfDataList,numeric-method}
\alias{plot,cfDataList,missing-method}
\alias{plot,cfOther,missing-method}
\title{Default Clifro Plotting}
\usage{
\S4method{plot}{cfDataList,numeric}(x, y, ...)

\S4method{plot}{cfDataList,missing}(x, y, ...)

\S4method{plot}{cfOther,missing}(x, y)
}
\arguments{
\item{x}{a \code{cfData} or \code{cfDataList} object.}

\item{y}{missing for \code{cfData} objects, or a number representing the
dataframe to plot if \code{x} is a \code{cfDataList} object.}

\item{...}{arguments passed onto the different plotting methods.

These methods are intended to simplify the data visualisation and exploration
of CliFlo data. The type of plot is determined by the type of the data output
from a \pkg{clifro} query. All of these methods plot individual plots for
each CliFlo station (if there is more than one in the query). If \code{x} is
a \code{cfDataList}, by default the first datatype will be plotted unless
\code{y} is supplied.

The following table links the datatypes to the corresponding plot methods:

\tabular{ll}{
\strong{Datatype} \tab \strong{Method}\cr
Wind \tab \code{\link{plot.cfWind}} for windrose, wind speed and direction
contour plots\cr
Rain \tab \code{\link{plot.cfRain}} for plotting rainfall (mm) through time\cr
Screen Obs \tab \code{\link{plot.cfScreenObs}} for time series plots of air,
wet bulb, and dew-point temperature plots\cr
Max/Min Temp \tab \code{\link{plot.cfTemp}} for maximum, minimum and
average temperature time series plots\cr
Earth Temp \tab \code{\link{plot.cfEarthTemp}} for earth temperature
time series plots\cr
Sunshine \tab \code{\link{plot.cfSunshine}} for accumulated, hourly or daily
sunshine, time series plots\cr
Pressure \tab \code{\link{plot.cfPressure}} for mean sea level atmospheric
pressure time series plots\cr
Other data \tab No default plot methods\cr
}}
}
\description{
Plot \pkg{clifro} data based on the datatype.
}
\seealso{
\code{\link{cf_query}} to retrieve the CliFlo data and create
 \code{cfData} objects.

  Refer to \code{\link[ggplot2]{theme}} for more possible arguments to pass
  to these methods.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findStations.R
\name{cf_find_station}
\alias{cf_find_station}
\title{Search for Clifro Stations}
\usage{
cf_find_station(
  ...,
  search = c("name", "region", "network", "latlong"),
  datatype,
  combine = c("all", "any"),
  status = c("open", "closed", "all")
)
}
\arguments{
\item{...}{arguments to pass into the search, these differ depending on
\code{search}.}

\item{search}{one of \code{name}, \code{network}, \code{region} or
\code{latlong} indicating the type of search to be conducted.}

\item{datatype}{\code{cfDatatype} object for when the search is based on
datatypes.}

\item{combine}{character string \code{"all"} or \code{"any"} indicating if the
stations contain all or any of the selected datatypes for when the search is
based on datatypes.}

\item{status}{character string indicating \code{"open"}, \code{"closed"} or
\code{"all"} stations be returned by the search.}
}
\value{
\code{cfStation} object
}
\description{
Search for \pkg{clifro} stations based on name, region, location or network
number, and return a \code{cfStation} object.
}
\details{
The \code{cf_find_station} function is a convenience function for finding
CliFlo stations in \R. It uses the CliFlo
\href{https://cliflo.niwa.co.nz/pls/niwp/wstn.get_stn_html}{Find Stations}
page to do the searching, and therefore means that the stations are not
stored within \pkg{clifro}.

If \code{datatype} is missing then the search is conducted
without any reference to datatypes. If it is supplied then the
search will only return stations that have any or all of the supplied
datatypes, depending on \code{combine}. The default behaviour is to search
for stations based on pattern matching the station name and return only the
open stations.

If the \code{latlong} search type is used the function expects named
arguments with names (partially) matching latitude,
longitude and radius. If the arguments are passed in without names they must
be in order of latitude, longitude and radius (see examples).
}
\note{
Since the searching is done by CliFlo there are obvious restrictions.
Unfortunately the pattern matching for station name does not provide
functionality for regular expressions, nor does it allow simultaneous
searches although \pkg{clifro} does provide some extra functionality, see
the 'OR query Search' example below.
}
\examples{
\dontrun{
# Station Name Search ------------------------------------------------------
# Return all open stations with 'island' in the name (pattern match search)
# Note this example uses all the defaults

island_st = cf_find_station("island")
island_st

# Region Search ------------------------------------------------------------
# Return all the closed stations from Queenstown (using partial matching)

queenstown.st = cf_find_station("queen", search = "region", status = "closed")
queenstown.st

# Long/Lat Search ----------------------------------------------------------
# Return all open stations within a 10km radius of the Beehive in Wellington
# From Wikipedia: latitude 41.2784 S, longitude 174.7767 E

beehive.st = cf_find_station(lat = -41.2784, long = 174.7767, rad = 10,
                             search = "latlong")
beehive.st

# Network ID Search --------------------------------------------------------
# Return all stations that share A42 in their network ID

A42.st = cf_find_station("A42", search = "network", status = "all")
A42.st

# Using Datatypes in the Search --------------------------------------------
# Is the Reefton EWS station open and does it collect daily rain and/or wind
# data?

# First, create the daily rain and wind datatypes
daily.dt = cf_datatype(c(2, 3), c(1, 1), list(4, 1), c(1, NA))
daily.dt

# Then combine into the search. This will only return stations where at least
# one datatype is available.
cf_find_station("reefton EWS", datatype = daily.dt)  # Yes

# OR Query Search ----------------------------------------------------------
# Return all stations sharing A42 in their network ID *or* all the open
# stations within 10km of the Beehive in Wellington (note this is not
# currently available as a single query in CliFlo).

cf_find_station("A42", search = "network", status = "all") +
cf_find_station(lat = -41.2784, long = 174.7767, rad = 10,
                search = "latlong")

# Note these are all ordered by open stations, then again by their end dates
}
}
\seealso{
\code{\link{cf_save_kml}} for saving the resulting stations as a KML
file, \code{\link{cf_station}} for creating \code{\link{cfStation}} objects
when the agent numbers are known, \code{vignette("choose-station")} for a
tutorial on finding \pkg{clifro} stations and \code{vignette("cfStation")}
for working with \code{\link{cfStation}} objects.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cfUser.R
\name{summary,cfUser-method}
\alias{summary,cfUser-method}
\title{Summarise User Information}
\usage{
\S4method{summary}{cfUser}(object)
}
\arguments{
\item{object}{an object of class \code{cfUser}.}
}
\description{
Show the subscription status for the \pkg{clifro} user
}
