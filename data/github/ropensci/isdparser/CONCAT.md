isdparser
=========



[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/isdparser)](https://cranchecks.info/pkgs/isdparser)
[![R-check](https://github.com/ropensci/isdparser/workflows/R-check/badge.svg)](https://github.com/ropensci/isdparser/actions?query=workflow%3AR-check)
[![codecov.io](https://codecov.io/github/ropensci/isdparser/coverage.svg?branch=master)](https://codecov.io/github/ropensci/isdparser?branch=master)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/isdparser?color=C9A115)](https://github.com/metacran/cranlogs.app)
[![cran version](http://www.r-pkg.org/badges/version/isdparser)](https://cran.r-project.org/package=isdparser)

Parse NOAA Integrated Surface Data Files

Documentation at https://docs.ropensci.org/isdparser/

## isdparser: Parse 'NOAA' Integrated Surface Data Files:
`isdparser` is a parser for 'NOAA' Integrated Surface Data ('ISD') files, described at <https://www.ncdc.noaa.gov/isd>. ISD includes numerous parameters such as wind speed and direction, wind gust, temperature, dew point, cloud data, sea level pressure, altimeter setting, station pressure, present weather, visibility, precipitation amounts for various time periods, snow depth, and various other elements as observed by each station. Data is stored as variable length ASCII character strings, with most fields optional. Included are tools for parsing entire files, or individual lines of data.

### Coverage
ISD is a global database, with data from approximately 35,000 stations worldwide, though the best spatial coverage is evident in North America, Europe, Australia, and parts of Asia. Coverage in the Northern Hemisphere is better than the Southern Hemisphere, and the overall period of record is currently 1901 to present. 


Code liberated from `rnoaa` to focus on ISD parsing since it's sorta
complicated. Has minimal dependencies, so you can parse your ISD/ISH
files without needing the deps that `rnoaa` needs. Will be used by
`rnoaa` once on CRAN.

Documentation at ftp://ftp.ncdc.noaa.gov/pub/data/noaa/ish-format-document.pdf

## API:

* `isd_parse()` - parse all lines in a file, with parallel option
* `isd_parse_line()` - parse a single line - you choose which lines to parse
and how to apply the function to your lines
* `isd_transform()` - transform ISD data variables
* `isd_parse_csv()` - parse csv format files

`isd_parse_csv()` parses NOAA ISD csv files, whereas `isd_parse()` and `isd_parse_line()`
both handle compressed files where each row of data is a string that needs to be parsed.

`isd_parse_csv()` is faster than `isd_parse()` because parsing each line takes some
time - although using `isd_parse(parallel = TRUE)` option gets closer to the speed of
`isd_parse_csv()`.

## Installation

CRAN stable version


```r
install.packages("isdparser")
```

Dev version


```r
remotes::install_github("ropensci/isdparser")
```


```r
library('isdparser')
```

## isd_parse_csv: parse a CSV file

Using a csv file included in the package:


```r
path <- system.file('extdata/00702699999.csv', package = "isdparser")
isd_parse_csv(path)
#> # A tibble: 6,843 x 68
#>    station date                source latitude longitude elevation name 
#>      <int> <dttm>               <int>    <dbl>     <dbl>     <dbl> <chr>
#>  1  7.03e8 2017-02-10 14:04:00      4        0         0      7026 WXPO…
#>  2  7.03e8 2017-02-10 14:14:00      4        0         0      7026 WXPO…
#>  3  7.03e8 2017-02-10 14:19:00      4        0         0      7026 WXPO…
#>  4  7.03e8 2017-02-10 14:24:00      4        0         0      7026 WXPO…
#>  5  7.03e8 2017-02-10 14:29:00      4        0         0      7026 WXPO…
#>  6  7.03e8 2017-02-10 14:34:00      4        0         0      7026 WXPO…
#>  7  7.03e8 2017-02-10 14:39:00      4        0         0      7026 WXPO…
#>  8  7.03e8 2017-02-10 14:44:00      4        0         0      7026 WXPO…
#>  9  7.03e8 2017-02-10 14:49:00      4        0         0      7026 WXPO…
#> 10  7.03e8 2017-02-10 14:54:00      4        0         0      7026 WXPO…
#> # … with 6,833 more rows, and 61 more variables: report_type <chr>,
#> #   call_sign <int>, quality_control <chr>, wnd <chr>, cig <chr>, vis <chr>,
#> #   tmp <chr>, dew <chr>, slp <chr>, wind_direction <chr>,
#> #   wind_direction_quality <chr>, wind_code <chr>, wind_speed <chr>,
#> #   wind_speed_quality <chr>, ceiling_height <chr>,
#> #   ceiling_height_quality <chr>, ceiling_height_determination <chr>,
#> #   ceiling_height_cavok <chr>, visibility_distance <chr>,
#> #   visibility_distance_quality <chr>, visibility_code <chr>,
#> #   visibility_code_quality <chr>, temperature <chr>,
#> #   temperature_quality <chr>, temperature_dewpoint <chr>,
#> #   temperature_dewpoint_quality <chr>, air_pressure <chr>,
#> #   air_pressure_quality <chr>, automated_atmospheric_condition_code <chr>,
#> #   quality_automated_atmospheric_condition_code <chr>, coverage_code <chr>,
#> #   coverage_quality_code <chr>, base_height_dimension <chr>,
#> #   base_height_quality_code <chr>, cloud_type_code <chr>,
#> #   cloud_type_quality_code <chr>, connective_cloud_attribute <chr>,
#> #   vertical_datum_attribute <chr>, base_height_upper_range_attribute <chr>,
#> #   base_height_lower_range_attribute <chr>, coverage <chr>,
#> #   opaque_coverage <chr>, coverage_quality <chr>, lowest_cover <chr>,
#> #   lowest_cover_quality <chr>, low_cloud_genus <chr>,
#> #   low_cloud_genus_quality <chr>, lowest_cloud_base_height <chr>,
#> #   lowest_cloud_base_height_quality <chr>, mid_cloud_genus <chr>,
#> #   mid_cloud_genus_quality <chr>, high_cloud_genus <chr>,
#> #   high_cloud_genus_quality <chr>, altimeter_setting_rate <chr>,
#> #   altimeter_quality_code <chr>, station_pressure_rate <chr>,
#> #   station_pressure_quality_code <chr>, speed_rate <chr>, quality_code <chr>,
#> #   rem <chr>, eqd <chr>
```

### isd_parse: parse a file with ASCII strings


```r
path <- system.file('extdata/024130-99999-2016.gz', package = "isdparser")
isd_parse(path)
#> # A tibble: 2,601 x 38
#>    total_chars usaf_station wban_station date  time  date_flag latitude
#>    <chr>       <chr>        <chr>        <chr> <chr> <chr>     <chr>   
#>  1 0054        024130       99999        2016… 0000  4         +60750  
#>  2 0054        024130       99999        2016… 0100  4         +60750  
#>  3 0054        024130       99999        2016… 0200  4         +60750  
#>  4 0054        024130       99999        2016… 0300  4         +60750  
#>  5 0054        024130       99999        2016… 0400  4         +60750  
#>  6 0039        024130       99999        2016… 0500  4         +60750  
#>  7 0054        024130       99999        2016… 0600  4         +60750  
#>  8 0039        024130       99999        2016… 0700  4         +60750  
#>  9 0054        024130       99999        2016… 0800  4         +60750  
#> 10 0054        024130       99999        2016… 0900  4         +60750  
#> # … with 2,591 more rows, and 31 more variables: longitude <chr>,
#> #   type_code <chr>, elevation <chr>, call_letter <chr>, quality <chr>,
#> #   wind_direction <chr>, wind_direction_quality <chr>, wind_code <chr>,
#> #   wind_speed <chr>, wind_speed_quality <chr>, ceiling_height <chr>,
#> #   ceiling_height_quality <chr>, ceiling_height_determination <chr>,
#> #   ceiling_height_cavok <chr>, visibility_distance <chr>,
#> #   visibility_distance_quality <chr>, visibility_code <chr>,
#> #   visibility_code_quality <chr>, temperature <chr>,
#> #   temperature_quality <chr>, temperature_dewpoint <chr>,
#> #   temperature_dewpoint_quality <chr>, air_pressure <chr>,
#> #   air_pressure_quality <chr>,
#> #   AW1_present_weather_observation_identifier <chr>,
#> #   AW1_automated_atmospheric_condition_code <chr>,
#> #   AW1_quality_automated_atmospheric_condition_code <chr>, REM_remarks <chr>,
#> #   REM_identifier <chr>, REM_length_quantity <chr>, REM_comment <chr>
```

process in parallel


```r
isd_parse(path, parallel = TRUE)
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/isdparser/issues).
* License: MIT
* Get citation information for `isdparser` in R doing `citation(package = 'isdparser')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![rofooter](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
isdparser 0.4.0
===============

### NEW FEATURES

* Gains new function `isd_parse_csv()` for working with ISD csv files rather than the compressed ASCII text files that `isd_parse()` parses. This new function introduces two new package imports: data.table and lubridate. The internals of the parsing tools for the various data types have changed to accommodate both csv and ASCII text files (#16)

isdparser 0.3.0
===============

### MINOR IMPROVEMENTS

* towards working on integrating metadata for each of the data fields, we've incorporated a data.frame of metadata into the package. see `?isd_metadata`. eventually we'd like to allow conversions and such on the data based on units (#12)

### BUG FIXES

* fix to internal parsing of string; first remove REM section from the ADD section so that codes in REM that happen to match those in ADD aren't detected spuriously (#15) thanks @ezwelty


isdparser 0.2.0
===============

### NEW FEATURES

* New function `isd_transform()` that transforms some variables in 
ISD output data - but not all. This functionality used to be done
by default in `isd_parse()` and `isd_parse_line()` - but data transforming 
yanked out as separate step done with `isd_transform()` (#11)
* Both `isd_parse()` and `isd_parse_line()` gain new parameter 
`additional` to optionally include additional and remarks
data sections in output (#10)

### BUG FIXES

* Fixed a bug in `isd()` in which a data section was not parsed 
correctly. (#9)


isdparser 0.1.0
===============

### NEW FEATURES

* Released to CRAN.
## Test environments

* local OS X install, R 3.6.2 Patched
* ubuntu 16.04 (on travis-ci), R 3.6.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

* Checked the 1 reverse dependency and there was no problem. 

--- 

This version adds a new function to parse a different data format from the same data source.

Thanks!
Scott Chamberlain
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

* Submit an issue on the [Issues page](https://github.com/ropensci/isdparser/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/isdparser.git`
* Make sure to track progress upstream (i.e., on our version of `isdparser` at `ropensci/isdparser`) by doing `git remote add upstream https://github.com/ropensci/isdparser.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Please do write a test(s) for your changes if they affect code and not just docs
* Push up to your account
* Submit a pull request to home base at `ropensci/isdparser`

### Also, check out our [discussion forum](https://discuss.ropensci.org)
<!-- Do not share screen shots of code. Share actual code in text format. -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/). If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
# Platform

|field    |value                                       |
|:--------|:-------------------------------------------|
|version  |R version 3.6.2 Patched (2020-01-25 r77715) |
|os       |macOS Mojave 10.14.6                        |
|system   |x86_64, darwin15.6.0                        |
|ui       |X11                                         |
|language |(EN)                                        |
|collate  |en_US.UTF-8                                 |
|ctype    |en_US.UTF-8                                 |
|tz       |US/Pacific                                  |
|date     |2020-01-31                                  |

# Dependencies

|package   |old   |new   |Δ  |
|:---------|:-----|:-----|:--|
|isdparser |0.3.0 |0.4.0 |*  |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*isdparser
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
[![cran checks](https://cranchecks.info/badges/worst/isdparser)](https://cranchecks.info/pkgs/isdparser)
[![R-check](https://github.com/ropensci/isdparser/workflows/R-check/badge.svg)](https://github.com/ropensci/isdparser/actions?query=workflow%3AR-check)
[![codecov.io](https://codecov.io/github/ropensci/isdparser/coverage.svg?branch=master)](https://codecov.io/github/ropensci/isdparser?branch=master)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/isdparser?color=C9A115)](https://github.com/metacran/cranlogs.app)
[![cran version](http://www.r-pkg.org/badges/version/isdparser)](https://cran.r-project.org/package=isdparser)

Parse NOAA Integrated Surface Data Files

Documentation at https://docs.ropensci.org/isdparser/

## isdparser: Parse 'NOAA' Integrated Surface Data Files:
`isdparser` is a parser for 'NOAA' Integrated Surface Data ('ISD') files, described at <https://www.ncdc.noaa.gov/isd>. ISD includes numerous parameters such as wind speed and direction, wind gust, temperature, dew point, cloud data, sea level pressure, altimeter setting, station pressure, present weather, visibility, precipitation amounts for various time periods, snow depth, and various other elements as observed by each station. Data is stored as variable length ASCII character strings, with most fields optional. Included are tools for parsing entire files, or individual lines of data.

### Coverage
ISD is a global database, with data from approximately 35,000 stations worldwide, though the best spatial coverage is evident in North America, Europe, Australia, and parts of Asia. Coverage in the Northern Hemisphere is better than the Southern Hemisphere, and the overall period of record is currently 1901 to present. 


Code liberated from `rnoaa` to focus on ISD parsing since it's sorta
complicated. Has minimal dependencies, so you can parse your ISD/ISH
files without needing the deps that `rnoaa` needs. Will be used by
`rnoaa` once on CRAN.

Documentation at ftp://ftp.ncdc.noaa.gov/pub/data/noaa/ish-format-document.pdf

## API:

* `isd_parse()` - parse all lines in a file, with parallel option
* `isd_parse_line()` - parse a single line - you choose which lines to parse
and how to apply the function to your lines
* `isd_transform()` - transform ISD data variables
* `isd_parse_csv()` - parse csv format files

`isd_parse_csv()` parses NOAA ISD csv files, whereas `isd_parse()` and `isd_parse_line()`
both handle compressed files where each row of data is a string that needs to be parsed.

`isd_parse_csv()` is faster than `isd_parse()` because parsing each line takes some
time - although using `isd_parse(parallel = TRUE)` option gets closer to the speed of
`isd_parse_csv()`.

## Installation

CRAN stable version

```{r eval=FALSE}
install.packages("isdparser")
```

Dev version

```{r eval=FALSE}
remotes::install_github("ropensci/isdparser")
```

```{r}
library('isdparser')
```

## isd_parse_csv: parse a CSV file

Using a csv file included in the package:

```{r}
path <- system.file('extdata/00702699999.csv', package = "isdparser")
isd_parse_csv(path)
```

### isd_parse: parse a file with ASCII strings

```{r}
path <- system.file('extdata/024130-99999-2016.gz', package = "isdparser")
isd_parse(path)
```

process in parallel

```{r eval=FALSE}
isd_parse(path, parallel = TRUE)
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/isdparser/issues).
* License: MIT
* Get citation information for `isdparser` in R doing `citation(package = 'isdparser')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![rofooter](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
---
title: "Introduction to the isdparser package"
author: "Scott Chamberlain"
date: "2020-01-31"
output:
  html_document:
    toc: true
    toc_float: true
vignette: >
  %\VignetteIndexEntry{Introduction to the isdparser package}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



`isdparser` is an parser for ISD/ISD NOAA files

Code liberated from `rnoaa` to focus on ISD parsing since it's sorta
complicated. Has minimal dependencies, so you can parse your ISD/ISH
files without needing the deps that `rnoaa` needs. Will be used by
`rnoaa` once on CRAN.

Documentation at ftp://ftp.ncdc.noaa.gov/pub/data/noaa/ish-format-document.pdf

Package API:

* `isd_parse()` - parse all lines in a file, with parallel option
* `isd_parse_line()` - parse a single line - you choose which lines to parse
and how to apply the function to your lines
* `isd_transform()` - transform ISD data variables
* `isd_parse_csv()` - parse csv format files

`isd_parse_csv()` parses NOAA ISD csv files, whereas `isd_parse()` and `isd_parse_line()`
both handle compressed files where each row of data is a string that needs to be parsed.

`isd_parse_csv()` is faster than `isd_parse()` because parsing each line takes some
time - although using `isd_parse(parallel = TRUE)` option gets closer to the speed of
`isd_parse_csv()`.


## Install

Stable from CRAN


```r
install.packages("isdparser")
```

Dev version


```r
remotes::install_github("ropensci/isdparser")
```


```r
library("isdparser")
```

## isd_parse_csv: parse a CSV file

Using a csv file included in the package:


```r
path <- system.file('extdata/00702699999.csv', package = "isdparser")
isd_parse_csv(path)
#> # A tibble: 6,843 x 68
#>    station date                source latitude longitude elevation name 
#>      <int> <dttm>               <int>    <dbl>     <dbl>     <dbl> <chr>
#>  1  7.03e8 2017-02-10 14:04:00      4        0         0      7026 WXPO…
#>  2  7.03e8 2017-02-10 14:14:00      4        0         0      7026 WXPO…
#>  3  7.03e8 2017-02-10 14:19:00      4        0         0      7026 WXPO…
#>  4  7.03e8 2017-02-10 14:24:00      4        0         0      7026 WXPO…
#>  5  7.03e8 2017-02-10 14:29:00      4        0         0      7026 WXPO…
#>  6  7.03e8 2017-02-10 14:34:00      4        0         0      7026 WXPO…
#>  7  7.03e8 2017-02-10 14:39:00      4        0         0      7026 WXPO…
#>  8  7.03e8 2017-02-10 14:44:00      4        0         0      7026 WXPO…
#>  9  7.03e8 2017-02-10 14:49:00      4        0         0      7026 WXPO…
#> 10  7.03e8 2017-02-10 14:54:00      4        0         0      7026 WXPO…
#> # … with 6,833 more rows, and 61 more variables: report_type <chr>,
#> #   call_sign <int>, quality_control <chr>, wnd <chr>, cig <chr>, vis <chr>,
#> #   tmp <chr>, dew <chr>, slp <chr>, wind_direction <chr>,
#> #   wind_direction_quality <chr>, wind_code <chr>, wind_speed <chr>,
#> #   wind_speed_quality <chr>, ceiling_height <chr>,
#> #   ceiling_height_quality <chr>, ceiling_height_determination <chr>,
#> #   ceiling_height_cavok <chr>, visibility_distance <chr>,
#> #   visibility_distance_quality <chr>, visibility_code <chr>,
#> #   visibility_code_quality <chr>, temperature <chr>,
#> #   temperature_quality <chr>, temperature_dewpoint <chr>,
#> #   temperature_dewpoint_quality <chr>, air_pressure <chr>,
#> #   air_pressure_quality <chr>, automated_atmospheric_condition_code <chr>,
#> #   quality_automated_atmospheric_condition_code <chr>, coverage_code <chr>,
#> #   coverage_quality_code <chr>, base_height_dimension <chr>,
#> #   base_height_quality_code <chr>, cloud_type_code <chr>,
#> #   cloud_type_quality_code <chr>, connective_cloud_attribute <chr>,
#> #   vertical_datum_attribute <chr>, base_height_upper_range_attribute <chr>,
#> #   base_height_lower_range_attribute <chr>, coverage <chr>,
#> #   opaque_coverage <chr>, coverage_quality <chr>, lowest_cover <chr>,
#> #   lowest_cover_quality <chr>, low_cloud_genus <chr>,
#> #   low_cloud_genus_quality <chr>, lowest_cloud_base_height <chr>,
#> #   lowest_cloud_base_height_quality <chr>, mid_cloud_genus <chr>,
#> #   mid_cloud_genus_quality <chr>, high_cloud_genus <chr>,
#> #   high_cloud_genus_quality <chr>, altimeter_setting_rate <chr>,
#> #   altimeter_quality_code <chr>, station_pressure_rate <chr>,
#> #   station_pressure_quality_code <chr>, speed_rate <chr>, quality_code <chr>,
#> #   rem <chr>, eqd <chr>
```

Download a file first:


```r
path <- file.path(tempdir(), "00702699999.csv")
x <- "https://www.ncei.noaa.gov/data/global-hourly/access/2017/00702699999.csv"
download.file(x, path)
isd_parse_csv(path)
#> # A tibble: 6,843 x 68
#>    station date                source latitude longitude elevation name 
#>      <int> <dttm>               <int>    <dbl>     <dbl>     <dbl> <chr>
#>  1  7.03e8 2017-02-10 14:04:00      4        0         0      7026 WXPO…
#>  2  7.03e8 2017-02-10 14:14:00      4        0         0      7026 WXPO…
#>  3  7.03e8 2017-02-10 14:19:00      4        0         0      7026 WXPO…
#>  4  7.03e8 2017-02-10 14:24:00      4        0         0      7026 WXPO…
#>  5  7.03e8 2017-02-10 14:29:00      4        0         0      7026 WXPO…
#>  6  7.03e8 2017-02-10 14:34:00      4        0         0      7026 WXPO…
#>  7  7.03e8 2017-02-10 14:39:00      4        0         0      7026 WXPO…
#>  8  7.03e8 2017-02-10 14:44:00      4        0         0      7026 WXPO…
#>  9  7.03e8 2017-02-10 14:49:00      4        0         0      7026 WXPO…
#> 10  7.03e8 2017-02-10 14:54:00      4        0         0      7026 WXPO…
#> # … with 6,833 more rows, and 61 more variables: report_type <chr>,
#> #   call_sign <int>, quality_control <chr>, wnd <chr>, cig <chr>, vis <chr>,
#> #   tmp <chr>, dew <chr>, slp <chr>, wind_direction <chr>,
#> #   wind_direction_quality <chr>, wind_code <chr>, wind_speed <chr>,
#> #   wind_speed_quality <chr>, ceiling_height <chr>,
#> #   ceiling_height_quality <chr>, ceiling_height_determination <chr>,
#> #   ceiling_height_cavok <chr>, visibility_distance <chr>,
#> #   visibility_distance_quality <chr>, visibility_code <chr>,
#> #   visibility_code_quality <chr>, temperature <chr>,
#> #   temperature_quality <chr>, temperature_dewpoint <chr>,
#> #   temperature_dewpoint_quality <chr>, air_pressure <chr>,
#> #   air_pressure_quality <chr>, automated_atmospheric_condition_code <chr>,
#> #   quality_automated_atmospheric_condition_code <chr>, coverage_code <chr>,
#> #   coverage_quality_code <chr>, base_height_dimension <chr>,
#> #   base_height_quality_code <chr>, cloud_type_code <chr>,
#> #   cloud_type_quality_code <chr>, connective_cloud_attribute <chr>,
#> #   vertical_datum_attribute <chr>, base_height_upper_range_attribute <chr>,
#> #   base_height_lower_range_attribute <chr>, coverage <chr>,
#> #   opaque_coverage <chr>, coverage_quality <chr>, lowest_cover <chr>,
#> #   lowest_cover_quality <chr>, low_cloud_genus <chr>,
#> #   low_cloud_genus_quality <chr>, lowest_cloud_base_height <chr>,
#> #   lowest_cloud_base_height_quality <chr>, mid_cloud_genus <chr>,
#> #   mid_cloud_genus_quality <chr>, high_cloud_genus <chr>,
#> #   high_cloud_genus_quality <chr>, altimeter_setting_rate <chr>,
#> #   altimeter_quality_code <chr>, station_pressure_rate <chr>,
#> #   station_pressure_quality_code <chr>, speed_rate <chr>, quality_code <chr>,
#> #   rem <chr>, eqd <chr>
```



## isd_parse_line: parse lines from an ASCII strings file


```r
path <- system.file('extdata/024130-99999-2016.gz', package = "isdparser")
lns <- readLines(path, encoding = "latin1")
isd_parse_line(lns[1])
#> # A tibble: 1 x 38
#>   total_chars usaf_station wban_station date  time  date_flag latitude longitude
#>   <chr>       <chr>        <chr>        <chr> <chr> <chr>     <chr>    <chr>    
#> 1 0054        024130       99999        2016… 0000  4         +60750   +012767  
#> # … with 30 more variables: type_code <chr>, elevation <chr>,
#> #   call_letter <chr>, quality <chr>, wind_direction <chr>,
#> #   wind_direction_quality <chr>, wind_code <chr>, wind_speed <chr>,
#> #   wind_speed_quality <chr>, ceiling_height <chr>,
#> #   ceiling_height_quality <chr>, ceiling_height_determination <chr>,
#> #   ceiling_height_cavok <chr>, visibility_distance <chr>,
#> #   visibility_distance_quality <chr>, visibility_code <chr>,
#> #   visibility_code_quality <chr>, temperature <chr>,
#> #   temperature_quality <chr>, temperature_dewpoint <chr>,
#> #   temperature_dewpoint_quality <chr>, air_pressure <chr>,
#> #   air_pressure_quality <chr>,
#> #   AW1_present_weather_observation_identifier <chr>,
#> #   AW1_automated_atmospheric_condition_code <chr>,
#> #   AW1_quality_automated_atmospheric_condition_code <chr>, REM_remarks <chr>,
#> #   REM_identifier <chr>, REM_length_quantity <chr>, REM_comment <chr>
```

Or, give back a list


```r
head(
  isd_parse_line(lns[1], as_data_frame = FALSE)
)
#> $total_chars
#> [1] "0054"
#> 
#> $usaf_station
#> [1] "024130"
#> 
#> $wban_station
#> [1] "99999"
#> 
#> $date
#> [1] "20160101"
#> 
#> $time
#> [1] "0000"
#> 
#> $date_flag
#> [1] "4"
```

Optionally don't include "Additional" and "Remarks" sections in parsed output.


```r
isd_parse_line(lns[1], additional = FALSE)
#> # A tibble: 1 x 31
#>   total_chars usaf_station wban_station date  time  date_flag latitude longitude
#>   <chr>       <chr>        <chr>        <chr> <chr> <chr>     <chr>    <chr>    
#> 1 0054        024130       99999        2016… 0000  4         +60750   +012767  
#> # … with 23 more variables: type_code <chr>, elevation <chr>,
#> #   call_letter <chr>, quality <chr>, wind_direction <chr>,
#> #   wind_direction_quality <chr>, wind_code <chr>, wind_speed <chr>,
#> #   wind_speed_quality <chr>, ceiling_height <chr>,
#> #   ceiling_height_quality <chr>, ceiling_height_determination <chr>,
#> #   ceiling_height_cavok <chr>, visibility_distance <chr>,
#> #   visibility_distance_quality <chr>, visibility_code <chr>,
#> #   visibility_code_quality <chr>, temperature <chr>,
#> #   temperature_quality <chr>, temperature_dewpoint <chr>,
#> #   temperature_dewpoint_quality <chr>, air_pressure <chr>,
#> #   air_pressure_quality <chr>
```


## isd_parse: parse an ASCII strings file

Downloading a new file


```r
path <- file.path(tempdir(), "007026-99999-2017.gz")
y <- "ftp://ftp.ncdc.noaa.gov/pub/data/noaa/2017/007026-99999-2017.gz"
download.file(y, path)
isd_parse(path)
#> # A tibble: 6,843 x 72
#>    total_chars usaf_station wban_station date  time  date_flag latitude
#>    <chr>       <chr>        <chr>        <chr> <chr> <chr>     <chr>   
#>  1 0157        007026       99999        2017… 1404  4         +00000  
#>  2 0157        007026       99999        2017… 1414  4         +00000  
#>  3 0157        007026       99999        2017… 1419  4         +00000  
#>  4 0157        007026       99999        2017… 1424  4         +00000  
#>  5 0157        007026       99999        2017… 1429  4         +00000  
#>  6 0144        007026       99999        2017… 1434  4         +00000  
#>  7 0157        007026       99999        2017… 1439  4         +00000  
#>  8 0157        007026       99999        2017… 1444  4         +00000  
#>  9 0172        007026       99999        2017… 1449  4         +00000  
#> 10 0157        007026       99999        2017… 1454  4         +00000  
#> # … with 6,833 more rows, and 65 more variables: longitude <chr>,
#> #   type_code <chr>, elevation <chr>, call_letter <chr>, quality <chr>,
#> #   wind_direction <chr>, wind_direction_quality <chr>, wind_code <chr>,
#> #   wind_speed <chr>, wind_speed_quality <chr>, ceiling_height <chr>,
#> #   ceiling_height_quality <chr>, ceiling_height_determination <chr>,
#> #   ceiling_height_cavok <chr>, visibility_distance <chr>,
#> #   visibility_distance_quality <chr>, visibility_code <chr>,
#> #   visibility_code_quality <chr>, temperature <chr>,
#> #   temperature_quality <chr>, temperature_dewpoint <chr>,
#> #   temperature_dewpoint_quality <chr>, air_pressure <chr>,
#> #   air_pressure_quality <chr>, GF1_sky_condition <chr>, GF1_coverage <chr>,
#> #   GF1_opaque_coverage <chr>, GF1_coverage_quality <chr>,
#> #   GF1_lowest_cover <chr>, GF1_lowest_cover_quality <chr>,
#> #   GF1_low_cloud_genus <chr>, GF1_low_cloud_genus_quality <chr>,
#> #   GF1_lowest_cloud_base_height <chr>,
#> #   GF1_lowest_cloud_base_height_quality <chr>, GF1_mid_cloud_genus <chr>,
#> #   GF1_mid_cloud_genus_quality <chr>, GF1_high_cloud_genus <chr>,
#> #   GF1_high_cloud_genus_quality <chr>, MA1_atmospheric_pressure <chr>,
#> #   MA1_altimeter_setting_rate <chr>, MA1_altimeter_quality_code <chr>,
#> #   MA1_station_pressure_rate <chr>, MA1_station_pressure_quality_code <chr>,
#> #   REM_remarks <chr>, REM_identifier <chr>, REM_length_quantity <chr>,
#> #   REM_comment <chr>, OC1_wind_gust_observation_identifier <chr>,
#> #   OC1_speed_rate <chr>, OC1_quality_code <chr>,
#> #   GA1_sky_cover_layer_identifier <chr>, GA1_coverage_code <chr>,
#> #   GA1_coverage_quality_code <chr>, GA1_base_height_dimension <chr>,
#> #   GA1_base_height_quality_code <chr>, GA1_cloud_type_code <chr>,
#> #   GA1_cloud_type_quality_code <chr>, GE1_sky_condition <chr>,
#> #   GE1_connective_cloud_attribute <chr>, GE1_vertical_datum_attribute <chr>,
#> #   GE1_base_height_upper_range_attribute <chr>,
#> #   GE1_base_height_lower_range_attribute <chr>,
#> #   AW1_present_weather_observation_identifier <chr>,
#> #   AW1_automated_atmospheric_condition_code <chr>,
#> #   AW1_quality_automated_atmospheric_condition_code <chr>
```

### Parallel


```r
isd_parse(path, parallel = TRUE)
```

### Progress

> note: Progress not printed if `parallel = TRUE`


```r
isd_parse(path, progress = TRUE)
#>
#>   |========================================================================================| 100%
#> # A tibble: 2,601 × 42
#>    total_chars usaf_station wban_station       date  time date_flag latitude longitude type_code
#>          <dbl>        <chr>        <chr>     <date> <chr>     <chr>    <dbl>     <dbl>     <chr>
#> 1           54       024130        99999 2016-01-01  0000         4    60.75    12.767     FM-12
#> 2           54       024130        99999 2016-01-01  0100         4    60.75    12.767     FM-12
#> 3           54       024130        99999 2016-01-01  0200         4    60.75    12.767     FM-12
#> 4           54       024130        99999 2016-01-01  0300         4    60.75    12.767     FM-12
#> 5           54       024130        99999 2016-01-01  0400         4    60.75    12.767     FM-12
#> 6           39       024130        99999 2016-01-01  0500         4    60.75    12.767     FM-12
#> 7           54       024130        99999 2016-01-01  0600         4    60.75    12.767     FM-12
#> 8           39       024130        99999 2016-01-01  0700         4    60.75    12.767     FM-12
#> 9           54       024130        99999 2016-01-01  0800         4    60.75    12.767     FM-12
#> 10          54       024130        99999 2016-01-01  0900         4    60.75    12.767     FM-12
#> # ... with 2,591 more rows, and 33 more variables: elevation <dbl>, call_letter <chr>, quality <chr>,
#> #   wind_direction <dbl>, wind_direction_quality <chr>, wind_code <chr>, wind_speed <dbl>,
#> #   wind_speed_quality <chr>, ceiling_height <chr>, ceiling_height_quality <chr>,
#> #   ceiling_height_determination <chr>, ceiling_height_cavok <chr>, visibility_distance <chr>,
#> #   visibility_distance_quality <chr>, visibility_code <chr>, visibility_code_quality <chr>,
#> #   temperature <dbl>, temperature_quality <chr>, temperature_dewpoint <dbl>,
#> #   temperature_dewpoint_quality <chr>, air_pressure <dbl>, air_pressure_quality <chr>,
#> #   AW1_present_weather_observation_identifier <chr>, AW1_automated_atmospheric_condition_code <chr>,
#> #   AW1_quality_automated_atmospheric_condition_code <chr>, N03_original_observation <chr>,
#> #   N03_original_value_text <chr>, N03_units_code <chr>, N03_parameter_code <chr>, REM_remarks <chr>,
#> #   REM_identifier <chr>, REM_length_quantity <chr>, REM_comment <chr>
```

### Additional data

Optionally don't include "Additional" and "Remarks" sections in parsed output.


```r
isd_parse(path, additional = FALSE)
#> # A tibble: 6,843 x 31
#>    total_chars usaf_station wban_station date  time  date_flag latitude
#>    <chr>       <chr>        <chr>        <chr> <chr> <chr>     <chr>   
#>  1 0157        007026       99999        2017… 1404  4         +00000  
#>  2 0157        007026       99999        2017… 1414  4         +00000  
#>  3 0157        007026       99999        2017… 1419  4         +00000  
#>  4 0157        007026       99999        2017… 1424  4         +00000  
#>  5 0157        007026       99999        2017… 1429  4         +00000  
#>  6 0144        007026       99999        2017… 1434  4         +00000  
#>  7 0157        007026       99999        2017… 1439  4         +00000  
#>  8 0157        007026       99999        2017… 1444  4         +00000  
#>  9 0172        007026       99999        2017… 1449  4         +00000  
#> 10 0157        007026       99999        2017… 1454  4         +00000  
#> # … with 6,833 more rows, and 24 more variables: longitude <chr>,
#> #   type_code <chr>, elevation <chr>, call_letter <chr>, quality <chr>,
#> #   wind_direction <chr>, wind_direction_quality <chr>, wind_code <chr>,
#> #   wind_speed <chr>, wind_speed_quality <chr>, ceiling_height <chr>,
#> #   ceiling_height_quality <chr>, ceiling_height_determination <chr>,
#> #   ceiling_height_cavok <chr>, visibility_distance <chr>,
#> #   visibility_distance_quality <chr>, visibility_code <chr>,
#> #   visibility_code_quality <chr>, temperature <chr>,
#> #   temperature_quality <chr>, temperature_dewpoint <chr>,
#> #   temperature_dewpoint_quality <chr>, air_pressure <chr>,
#> #   air_pressure_quality <chr>
```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/isd_transform.R
\name{isd_transform}
\alias{isd_transform}
\title{Transform ISD data variables}
\usage{
isd_transform(x)
}
\arguments{
\item{x}{(data.frame/tbl_df) data.frame/tbl from \code{\link{isd_parse}} or
data.frame/tbl or list from \code{\link{isd_parse_line}}}
}
\value{
A tibble (data.frame) or list
}
\description{
Transform ISD data variables
}
\details{
This function helps you clean your ISD data. \code{\link{isd_parse}}
and \code{\link{isd_parse_line}} give back data without modifying the
data. However, you'll likely want to transform some of the variables,
in terms of the variable class (character to numeric), accounting for the
scaling factor (variable X may need to be multiplied by 1000 according
to the ISD docs), and missing values (unfortunately, missing value
standards vary across ISD data).
}
\section{operations performed}{

\itemize{
 \item scale latitude by factor of 1000
 \item scale longitude by factor of 1000
 \item scale elevation by factor of 10
 \item scale wind speed by factor of 10
 \item scale temperature by factor of 10
 \item scale temperature dewpoint by factor of 10
 \item scale air pressure by factor of 10
 \item scale precipitation by factor of 10
 \item convert date to a Date class with \code{\link{as.Date}}
 \item change wind direction to numeric
 \item change total characters to numeric
}
}

\examples{
path <- system.file('extdata/104270-99999-1928.gz', package = "isdparser")
(res <- isd_parse(path))
isd_transform(res)

lns <- readLines(path, encoding = "latin1")
# data.frame
(res <- isd_parse_line(lns[1]))
isd_transform(res)
# list
(res <- isd_parse_line(lns[1], as_data_frame = FALSE))
isd_transform(res)
}
\seealso{
\code{\link{isd_parse}}, \code{\link{isd_parse_line}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/isdparser-package.R
\docType{data}
\name{isd_metadata}
\alias{isd_metadata}
\title{NOAA ISD metadata data.frame}
\format{
A data frame with 643 rows and 19 columns
}
\description{
This data.frame includes metadata describing all the data provided in ISD
data files. And is used for transforming and scaling variables.
}
\details{
Original csv data is in inst/extdata/isd_metadata.csv, collected from


The data.frame has the following columns:

\itemize{
 \item pos - (chr) position, if any
 \item category - (chr) category, one of additional-data section,
 control-data section, element quality data section, mandatory-data section,
 original observation data section, or remarks data section
 \item sub_category - (chr) sub category label, one of climate reference
 network unique data, cloud and solar data, ground surface data, hail data,
 marine data, network metadata, precipitation-data, pressure data,
 runway visual range data, sea surface temperature, soil temperature data,
 temperature data, weather occurrence data, weather-occurrence-data,
 or wind data
 \item abbrev - (chr) abbreviation, if any, NA for control and mandatory
 sections
 \item label - (chr) label, a top level label for the data, usually the
 same as the abbreviation
 \item sub_label - (chr) sub label, a more detailed label about the
 variable
 \item field_length - (int) field length, number of characters
 \item min - (chr) minimum value, if applicable, original
 \item min_numeric - (int) minimum value, if applicable, integer
 \item max - (chr) maximum value, if applicable, original
 \item max_numeric - (chr) maximum value, if applicable, integer
 \item units - (chr) units, if applicable
 \item scaling_factor - (chr) scaling factor, original
 \item scaling_factor_numeric - (int) scaling factor, integer, one of
 1, 10, 100, 1000, or NA
 \item missing - (chr) value used to indicate missing data, original
 \item missing_numeric - (int) value used to indicate missing data, integer,
 one of 9, 99, 999, 9999, 99999, 999999, or NA
 \item description - (chr) short description of variable
 \item dom - (chr) long description of variable with categories
 \item dom_parsed_json - (list) NA if no categories, or a named list with
 category labels and their values
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/isd_parse_csv.R
\name{isd_parse_csv}
\alias{isd_parse_csv}
\title{Parse NOAA ISD/ISH csv data files}
\usage{
isd_parse_csv(path)
}
\arguments{
\item{path}{(character) file path. required}
}
\value{
A tibble (data.frame)
}
\description{
Parse NOAA ISD/ISH csv data files
}
\details{
Note that the `rem` (remarks) and `eqd` columns are
not parsed, just as with [isd_parse()].
}
\section{Column information}{


- USAF MASTER and NCEI WBAN station identifiers are combined into an 11
character code with the column `station`
- Date and Time have been combined to the column `date`
- Call letter is synonymous with `call_sign` column
- WIND-OBSERVATION is abbreviated as column `wnd`
- SKY-CONDITION-OBSERVATION is abbreviated as column `cig`
- VISIBILITY-OBSERVATION is abbreviated as column `vis`
- AIR-TEMPERATURE-OBSERVATION air temperature is abbreviated as the column
header `tmp`
- AIR-TEMPERATURE-OBSERVATION dew point is abbreviated as the column
`dew`
- AIR-PRESSURE-OBSERVATION sea level pressure is abbreviated as the column
`slp`
}

\examples{
path <- system.file('extdata/00702699999.csv', package = "isdparser")
(res <- isd_parse_csv(path))

# isd_parse_csv compared to isd_parse
if (interactive()) {
x="https://www.ncei.noaa.gov/data/global-hourly/access/2017/00702699999.csv"
download.file(x, (f_csv=file.path(tempdir(), "00702699999.csv")))
y="ftp://ftp.ncdc.noaa.gov/pub/data/noaa/2017/007026-99999-2017.gz"
download.file(y, (f_gz=file.path(tempdir(), "007026-99999-2017.gz")))
from_csv <- isd_parse_csv(f_csv)
from_gz <- isd_parse(f_gz, parallel = TRUE)

x="https://www.ncei.noaa.gov/data/global-hourly/access/1913/02982099999.csv"
download.file(x, (f=file.path(tempdir(), "02982099999.csv")))
isd_parse_csv(f)

x="https://www.ncei.noaa.gov/data/global-hourly/access/1923/02970099999.csv"
download.file(x, (f=file.path(tempdir(), "02970099999.csv")))
isd_parse_csv(f)

x="https://www.ncei.noaa.gov/data/global-hourly/access/1945/04390099999.csv"
download.file(x, (f=file.path(tempdir(), "04390099999.csv")))
isd_parse_csv(f)

x="https://www.ncei.noaa.gov/data/global-hourly/access/1976/02836099999.csv"
download.file(x, (f=file.path(tempdir(), "02836099999.csv")))
isd_parse_csv(f)
}
}
\references{
https://www.ncei.noaa.gov/data/global-hourly/access/
https://www.ncei.noaa.gov/data/global-hourly/doc/CSV_HELP.pdf
https://www.ncei.noaa.gov/data/global-hourly/doc/isd-format-document.pdf
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/isd_parser_line.R
\name{isd_parse_line}
\alias{isd_parse_line}
\title{Parse NOAA ISD/ISH data files - line by line}
\usage{
isd_parse_line(x, additional = TRUE, as_data_frame = TRUE)
}
\arguments{
\item{x}{(character) a single ISD line}

\item{additional}{(logical) include additional and remarks data sections
in output. Default: \code{TRUE}}

\item{as_data_frame}{(logical) output a tibble. Default: \code{FALSE}}
}
\value{
A tibble (data.frame)
}
\description{
Parse NOAA ISD/ISH data files - line by line
}
\examples{
path <- system.file('extdata/024130-99999-2016.gz', package = "isdparser")
lns <- readLines(path, encoding = "latin1")
isd_parse_line(lns[1])
isd_parse_line(lns[1], FALSE)

res <- lapply(lns[1:1000], isd_parse_line)
library("data.table")
library("tibble")
as_tibble(
 rbindlist(res, use.names = TRUE, fill = TRUE)
)

# only control + mandatory sections
isd_parse_line(lns[10], additional = FALSE)
isd_parse_line(lns[10], additional = TRUE)
}
\references{
ftp://ftp.ncdc.noaa.gov/pub/data/noaa
}
\seealso{
\code{\link{isd_parse}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/isdparser-package.R
\docType{package}
\name{isdparser-package}
\alias{isdparser-package}
\alias{isdparser}
\title{Parse NOAA ISD Files}
\section{Data format}{

Each record (data.frame row or individual list element) you get via
\code{isd_parse} or \code{isd_parse_line} has all data combined.
Control data fields are first, then mandatory fields, then additional data
fields and remarks. Control and mandatory fields have column names
describing what they are, while additional data fields have a length
three character prefix (e.g., AA1) linking the fields to the documentation
for the \strong{Additional Data Section} at
ftp://ftp.ncdc.noaa.gov/pub/data/noaa/ish-format-document.pdf
}

\section{Data size}{

Each line of an ISD data file has maximum of 2,844 characters.
}

\section{Control Data}{

The beginning of each record provides information about the report
including date, time, and station location information. Data fields
will be in positions identified in the applicable data definition.
The control data section is fixed length and is 60 characters long.
}

\section{Mandatory data}{

Each line of an ISD data file starts with mandatory data section.
The mandatory data section contains meteorological information on the
basic elements such as winds, visibility, and temperature. These are the
most commonly reported parameters and are available most of the time.
The mandatory data section is fixed length and is 45 characters long.
}

\section{Additional data}{

Each line of an ISD data file has an optional additional data
section, which follows the mandatory data section. These additional data
contain information of significance and/or which are received with
varying degrees of frequency. Identifiers are used to note when data
are present in the record. If all data fields in a group are missing,
the entire group is usually not reported. If no groups are reported
the section will be omitted. The additional data section is variable
in length with a minimum of 0 characters and a maximum of 637
(634 characters plus a 3 character section identifier) characters.
}

\section{Remarks data}{

The numeric and character (plain language) remarks are provided if they
exist. The data will vary in length and are identified in the applicable
data definition. The remarks section has a maximum length of 515
(512 characters plus a 3 character section identifier) characters.
}

\section{Missing values}{

Missing values for any non-signed item are filled (i.e., 999). Missing
values for any signed item are positive filled (i.e., +99999).
}

\section{Longitude and Latitude Coordinates}{

Longitudes will be reported with negative values representing longitudes
west of 0 degrees, and latitudes will be negative south of the equator.
Although the data field allows for values to a thousandth of a degree,
the values are often only computed to the hundredth of a degree with
a 0 entered in the thousandth position.
}

\author{
Scott Chamberlain \email{myrmecocystus@gmail.com}
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/isd_parser.R
\name{isd_parse}
\alias{isd_parse}
\title{Parse NOAA ISD/ISH data files}
\usage{
isd_parse(
  path,
  additional = TRUE,
  parallel = FALSE,
  cores = getOption("cl.cores", 2),
  progress = FALSE
)
}
\arguments{
\item{path}{(character) file path. required}

\item{additional}{(logical) include additional and remarks data sections
in output. Default: \code{TRUE}}

\item{parallel}{(logical). do processing in parallel. Default: \code{FALSE}}

\item{cores}{(integer) number of cores to use: Default: 2. We look in
your option "cl.cores", but use default value if not found.}

\item{progress}{(logical) print progress - ignored if \code{parallel=TRUE}.
The default is \code{FALSE} because printing progress adds a small bit of
time, so if processing time is important, then keep as \code{FALSE}}
}
\value{
A tibble (data.frame)
}
\description{
Parse NOAA ISD/ISH data files
}
\examples{
path <- system.file('extdata/104270-99999-1928.gz', package = "isdparser")

(res <- isd_parse(path))

# with progress
(res2 <- isd_parse(path, progress = TRUE))

# only control + mandatory sections
(res <- isd_parse(path, additional = FALSE))

\dontrun{
# in parallel
(out <- isd_parse(path, parallel = TRUE))
}
}
\references{
ftp://ftp.ncdc.noaa.gov/pub/data/noaa
}
\seealso{
\code{\link{isd_parse_line}}
}
