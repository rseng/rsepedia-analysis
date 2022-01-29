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


