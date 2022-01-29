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

*Wow, no problems at all. :)**Wow, no problems at all. :)*